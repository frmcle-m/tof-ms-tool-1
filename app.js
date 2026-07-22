/* ============================================================
   TOF-MS データ解析ツール — 実装コード
   (このファイルの中身は index.html の画面には表示されません。
    「コードを編集する」ページ（editor.html）から確認・編集できます)
   ============================================================ */

/* ============================================================
   1. 純粋な計算・パース関数（DOMに依存しない）
   ============================================================ */

// .887 測定条件ファイルのパース
function parseConditions(text) {
  const cond = {};
  const lines = text.split(/\r?\n/);
  for (let raw of lines) {
    const line = raw.trim();
    if (!line || line.startsWith(';')) continue;
    if (line.startsWith('REPORT-FILE')) {
      cond['measurement_time'] = line.replace('REPORT-FILE from ', '').split(' written')[0].trim();
      continue;
    }
    if (line.includes('=')) {
      const idx = line.indexOf('=');
      cond[line.slice(0, idx).trim()] = line.slice(idx + 1).trim();
      continue;
    }
    if (line.includes(':')) {
      const idx = line.indexOf(':');
      cond[line.slice(0, idx).trim()] = line.slice(idx + 1).trim();
      continue;
    }
    const parts = line.split(/\s+/);
    if (parts.length === 2) {
      cond[parts[0]] = parts[1];
    }
  }
  return cond;
}

// csv（タブ区切り: time_ns, counts）のパース
function parseCSV(text) {
  const lines = text.split(/\r?\n/);
  const time_ns = [];
  const counts = [];
  for (const line of lines) {
    if (!line) continue;
    let parts = line.split('\t');
    if (parts.length < 2) parts = line.trim().split(/\s+/);
    if (parts.length < 2) continue;
    const t = parseFloat(parts[0]);
    const c = parseFloat(parts[1]);
    if (Number.isNaN(t) || Number.isNaN(c)) continue;
    time_ns.push(t);
    counts.push(c);
  }
  return { time_ns: Float64Array.from(time_ns), counts: Float64Array.from(counts) };
}

// 測定条件の caloff を加えた時間軸 (ns) を作る。calfact は二重変換になるため使わない。
function makeTimeAxis(csvData, cond) {
  const caloff = parseFloat(cond && cond.caloff) || 0;
  return Float64Array.from(csvData.time_ns, v => v + caloff);
}

// numpy.convolve(..., mode='same') と同じ挙動の移動平均（端はゼロパディング）
function smooth(arr, window) {
  const w = Math.max(1, Math.floor(window) || 1);
  const n = arr.length;
  if (w <= 1) return Float64Array.from(arr);
  const out = new Float64Array(n);
  const prefix = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) prefix[i + 1] = prefix[i] + arr[i];
  const start = Math.floor((w - 1) / 2);
  for (let idx = 0; idx < n; idx++) {
    const fullIdx = idx + start;
    const xStart = fullIdx - (w - 1);
    const xEnd = fullIdx;
    const loClip = Math.max(xStart, 0);
    const hiClip = Math.min(xEnd, n - 1);
    let sum = 0;
    if (hiClip >= loClip) sum = prefix[hiClip + 1] - prefix[loClip];
    out[idx] = sum / w;
  }
  return out;
}

// scipy.signal.find_peaks(height=..., distance=...) 相当のピーク検出
function findPeaks(y, height, distance) {
  const dist = Math.max(1, Math.round(distance) || 1);
  const n = y.length;
  const candidates = [];
  for (let i = 1; i < n - 1; i++) {
    if (y[i] > y[i - 1] && y[i] > y[i + 1] && y[i] >= height) candidates.push(i);
  }
  candidates.sort((a, b) => y[b] - y[a]); // 高い順
  const kept = [];
  for (const idx of candidates) {
    let ok = true;
    for (const k of kept) {
      if (Math.abs(k - idx) < dist) { ok = false; break; }
    }
    if (ok) kept.push(idx);
  }
  kept.sort((a, b) => a - b);
  return kept;
}

// 最小二乗直線回帰: y = slope*x + intercept
function linreg(x, y) {
  const n = x.length;
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i];
  }
  const denom = n * sxx - sx * sx;
  const slope = denom !== 0 ? (n * sxy - sx * sy) / denom : 0;
  const intercept = (sy - slope * sx) / n;
  const meanY = sy / n;
  let ssTot = 0, ssRes = 0;
  for (let i = 0; i < n; i++) {
    const pred = slope * x[i] + intercept;
    ssRes += (y[i] - pred) ** 2;
    ssTot += (y[i] - meanY) ** 2;
  }
  const r2 = ssTot !== 0 ? 1 - ssRes / ssTot : 1;
  return { slope, intercept, r2 };
}

// 最適な開始質量(n)を n=1..50 で探索する質量校正
function calibrate(peakTimesS, massStep) {
  let best = { r2: -Infinity, slope: 0, intercept: 0, startN: 1 };
  for (let n = 1; n <= 50; n++) {
    const sqrtM = peakTimesS.map((_, i) => Math.sqrt(massStep * (n + i)));
    const { slope, intercept, r2 } = linreg(sqrtM, peakTimesS);
    if (r2 > best.r2) best = { r2, slope, intercept, startN: n };
  }
  return best;
}

// 時間(s) -> 質量(u)。t_eff <= 0 の点は mask=0（無効）
function timeToMass(tS, a, b) {
  const n = tS.length;
  const mass = new Float64Array(n);
  const mask = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    const tEff = tS[i] - b;
    if (tEff > 0) {
      mass[i] = (tEff / a) ** 2;
      mask[i] = 1;
    }
  }
  return { mass, mask };
}

// np.interp 相当の線形補間（xOldは昇順であること）
function interpLinear(xNew, xOld, yOld) {
  const n = xOld.length;
  const m = xNew.length;
  const out = new Float64Array(m);
  if (n === 0) return out;
  for (let i = 0; i < m; i++) {
    const x = xNew[i];
    if (x <= xOld[0]) { out[i] = yOld[0]; continue; }
    if (x >= xOld[n - 1]) { out[i] = yOld[n - 1]; continue; }
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (xOld[mid] <= x) lo = mid; else hi = mid;
    }
    const x0 = xOld[lo], x1 = xOld[hi], y0 = yOld[lo], y1 = yOld[hi];
    const t = x1 !== x0 ? (x - x0) / (x1 - x0) : 0;
    out[i] = y0 + t * (y1 - y0);
  }
  return out;
}

// ファイル名(stem)先頭の日付らしき数字列を抽出する（例: "20260625_4" -> "20260625"）。
// マッチしない場合は「未分類」の1グループにまとめる。
function extractDateKey(stem) {
  const m = stem.match(/^(\d{4,8})[-_]/);
  return m ? m[1] : '(日付不明)';
}

// 自然順(数値を含む文字列の番号順)ソート比較関数
// 例: "20260625_2" は "20260625_10" より前になる（文字列順だと逆転してしまうため）
function naturalCompare(a, b) {
  const re = /(\d+)|(\D+)/g;
  const aParts = a.match(re) || [];
  const bParts = b.match(re) || [];
  const len = Math.max(aParts.length, bParts.length);
  for (let i = 0; i < len; i++) {
    const ap = aParts[i] ?? '';
    const bp = bParts[i] ?? '';
    const aIsNum = /^\d+$/.test(ap);
    const bIsNum = /^\d+$/.test(bp);
    if (aIsNum && bIsNum) {
      const an = parseInt(ap, 10), bn = parseInt(bp, 10);
      if (an !== bn) return an - bn;
      if (ap !== bp) return ap < bp ? -1 : 1;
    } else if (ap !== bp) {
      return ap < bp ? -1 : 1;
    }
  }
  return 0;
}

function arraysEqual(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

function maskRange(arr, start, end) {
  const n = arr.length;
  const mask = new Uint8Array(n);
  for (let i = 0; i < n; i++) if (arr[i] >= start && arr[i] <= end) mask[i] = 1;
  return mask;
}

function sumMasked(arr, mask) {
  let s = 0;
  for (let i = 0; i < arr.length; i++) if (mask[i]) s += arr[i];
  return s;
}

function filterByMaskPair(xs, ys, mask) {
  const fx = [], fy = [];
  for (let i = 0; i < xs.length; i++) if (mask[i]) { fx.push(xs[i]); fy.push(ys[i]); }
  return { fx, fy };
}

const KEY_LABELS = {
  measurement_time: '測定時間',
  REALTIME: 'Real time (s)',
  TOTALSUM: 'Total counts',
  SWEEPS: 'Sweeps',
  swpreset: 'Sweep preset',
  range: 'Range (ch)',
  calunit: 'Cal. unit',
  calfact: 'Cal. factor',
  caloff: 'Cal. offset',
  syncout: 'Sync out (Hz)',
};
// calfact は整数値のため、小数点以下は切り捨てて表示する
function formatCalfact(v) {
  if (v === undefined || v === null || v === '') return v;
  const n = parseFloat(v);
  if (!Number.isFinite(n)) return v;
  return String(Math.trunc(n));
}
function makeConditionText(cond) {
  const lines = [];
  for (const key in KEY_LABELS) {
    if (!(key in cond)) continue;
    const value = key === 'calfact' ? formatCalfact(cond[key]) : cond[key];
    lines.push(`${KEY_LABELS[key]}: ${value}`);
  }
  return lines.join('\n');
}

/* ============================================================
   以下はブラウザ環境（DOM）でのみ実行する
   ============================================================ */
if (typeof document !== 'undefined') {

  /* ============================================================
     0. グローバル状態
     files: stem -> { name, csvFile, condFile, csvData, condData, label, labelType,
                       labelValue, dateKey, useInPlot, bgOverrideStem }
     knownDates:    これまでに一度でも読み込んだことのある日付キーの集合（サイドバーの日付一覧＝履歴）
     selectedDates: サイドバーの日付一覧で選ばれ、「2. 各ファイルの一覧・ラベル」に表示されている日付キーの集合（複数可）
     bgStem:        「共通BGとして使用」で指定された、共通のバックグラウンドファイル（ツール全体で1つ）
     checklistSelection: 複数ファイル重ね書き・積算タブのチェックリストの選択状態（stem -> true/false）
     fields:        表示設定・軸範囲・質量校正の入力値など（id -> 値）を記憶する
     ※「3. グラフを描画する」は日付を問わず全ファイルが対象のため、これらの設定はすべてツール全体で1つだけ持つ。
     ============================================================ */
  const state = {
    files: {},
    knownDates: new Set(),
    selectedDates: new Set(),
    bgStem: null,
    checklistSelection: { multiFileChecks: {}, intFileChecks: {} },
    fields: {},
  };

  function registerDateKey(key) {
    if (!state.knownDates.has(key)) {
      state.knownDates.add(key);
      state.selectedDates.add(key); // 新しく出現した日付は自動的に一覧表示の対象にする
    }
  }

  /* ---------- 1. ファイル読込ヘルパー ---------- */

  function extOf(name) {
    const m = name.match(/\.([^.]+)$/);
    return m ? m[1].toLowerCase() : '';
  }
  function stemOf(name) {
    return name.replace(/\.[^.]+$/, '');
  }
  function readFileAsText(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve(reader.result);
      reader.onerror = reject;
      reader.readAsText(file);
    });
  }
  async function ensureCsvParsed(stem) {
    const entry = state.files[stem];
    if (!entry || !entry.csvFile) return null;
    if (!entry.csvData) {
      const text = await readFileAsText(entry.csvFile);
      entry.csvData = parseCSV(text);
    }
    return entry.csvData;
  }

  /* ---------- 2. ローカル永続化 (IndexedDB) ----------
     読み込んだファイル（測定データ含む）と、表示設定・ラベル・BG設定をブラウザの
     IndexedDBに保存し、次回このHTMLを開いたときに自動で復元する。データは常にこの
     ブラウザ内にだけ保存され、外部サーバーには一切送信されない。 */
  const DB_NAME = 'tof_tool_db';
  const DB_VERSION = 4; // v4: ゴミ箱（trashストア）に対応
  const STORE_FILES = 'files';
  const STORE_META = 'meta';
  const STORE_TRASH = 'trash'; // ゴミ箱ページ（trash.html）からも同じ名前でアクセスする

  function openDB() {
    return new Promise((resolve) => {
      if (typeof indexedDB === 'undefined') { resolve(null); return; } // 非対応環境では永続化なしで動作継続
      let req;
      try {
        req = indexedDB.open(DB_NAME, DB_VERSION);
      } catch (e) {
        resolve(null);
        return;
      }
      req.onupgradeneeded = () => {
        const db = req.result;
        if (!db.objectStoreNames.contains(STORE_FILES)) db.createObjectStore(STORE_FILES, { keyPath: 'stem' });
        if (!db.objectStoreNames.contains(STORE_META)) db.createObjectStore(STORE_META, { keyPath: 'key' });
        if (!db.objectStoreNames.contains(STORE_TRASH)) db.createObjectStore(STORE_TRASH, { keyPath: 'stem' });
      };
      req.onsuccess = () => resolve(req.result);
      req.onerror = () => resolve(null);
    });
  }

  const dbPromise = openDB();
  if (typeof navigator !== 'undefined' && navigator.storage && navigator.storage.persist) {
    navigator.storage.persist().catch(() => {}); // ブラウザに保存領域の自動削除を避けるよう依頼する（対応環境のみ）
  }

  async function dbPut(storeName, value) {
    const db = await dbPromise;
    if (!db) return;
    return new Promise((resolve) => {
      const tx = db.transaction(storeName, 'readwrite');
      tx.objectStore(storeName).put(value);
      tx.oncomplete = () => resolve();
      tx.onerror = () => resolve();
    });
  }

  async function dbGetAll(storeName) {
    const db = await dbPromise;
    if (!db) return [];
    return new Promise((resolve) => {
      const tx = db.transaction(storeName, 'readonly');
      const req = tx.objectStore(storeName).getAll();
      req.onsuccess = () => resolve(req.result || []);
      req.onerror = () => resolve([]);
    });
  }

  async function dbDeleteKey(storeName, key) {
    const db = await dbPromise;
    if (!db) return;
    return new Promise((resolve) => {
      const tx = db.transaction(storeName, 'readwrite');
      tx.objectStore(storeName).delete(key);
      tx.oncomplete = () => resolve();
      tx.onerror = () => resolve();
    });
  }

  async function dbClear() {
    const db = await dbPromise;
    if (!db) return;
    return new Promise((resolve) => {
      const tx = db.transaction([STORE_FILES, STORE_META, STORE_TRASH], 'readwrite');
      tx.objectStore(STORE_FILES).clear();
      tx.objectStore(STORE_META).clear();
      tx.objectStore(STORE_TRASH).clear();
      tx.oncomplete = () => resolve();
      tx.onerror = () => resolve();
    });
  }

  // ファイルの状態を、IndexedDBに保存できる形（プレーンオブジェクト）に変換する
  function fileEntryToRecord(stem) {
    const e = state.files[stem];
    if (!e) return null;
    return {
      stem,
      csvData: e.csvData ? { time_ns: Array.from(e.csvData.time_ns), counts: Array.from(e.csvData.counts) } : null,
      condData: e.condData || {},
      label: e.label, labelType: e.labelType, labelValue: e.labelValue,
      useInPlot: e.useInPlot, dateKey: e.dateKey,
      hasCsv: !!e.csvFile, hasCond: !!e.condFile,
      bgOverrideStem: e.bgOverrideStem || null,
    };
  }

  // 1ファイル分（測定データ・条件・ラベル・個別BG指定等）をIndexedDBに保存する
  function persistFileEntry(stem) {
    const record = fileEntryToRecord(stem);
    if (!record) return;
    dbPut(STORE_FILES, record);
  }

  // ファイル単位ではない設定（表示設定・BG・日付一覧の選択状態など）をIndexedDBに保存する
  // （入力のたびに呼ばれるため、軽く間引く）
  let persistMetaTimer = null;
  function persistMetaNow() {
    dbPut(STORE_META, {
      key: 'settings',
      knownDates: [...state.knownDates],
      selectedDates: [...state.selectedDates],
      bgStem: state.bgStem,
      checklistSelection: state.checklistSelection,
      fields: state.fields,
    });
  }
  function schedulePersistMeta() {
    if (persistMetaTimer) clearTimeout(persistMetaTimer);
    persistMetaTimer = setTimeout(persistMetaNow, 300);
  }

  // 起動時にIndexedDBから前回のファイル・設定を読み戻す
  async function restoreFromDB() {
    const records = await dbGetAll(STORE_FILES);
    for (const r of records) {
      state.files[r.stem] = {
        name: r.stem,
        csvFile: r.hasCsv ? true : null,   // 実ファイルの代わりに真偽値のプレースホルダを置く（csvDataは既に復元済みのため不要）
        condFile: r.hasCond ? true : null,
        csvData: r.csvData ? { time_ns: Float64Array.from(r.csvData.time_ns), counts: Float64Array.from(r.csvData.counts) } : null,
        condData: r.condData || {},
        label: r.label || '', labelType: r.labelType || null,
        labelValue: r.labelValue != null ? r.labelValue : null,
        useInPlot: r.useInPlot !== false, dateKey: r.dateKey || extractDateKey(r.stem),
        bgOverrideStem: r.bgOverrideStem || null,
      };
      registerDateKey(state.files[r.stem].dateKey);
    }
    const metaRecords = await dbGetAll(STORE_META);
    const meta = metaRecords.find(m => m.key === 'settings');
    if (meta) {
      if (meta.knownDates) state.knownDates = new Set(meta.knownDates);
      if (meta.selectedDates) state.selectedDates = new Set(meta.selectedDates);
      if (meta.fields) state.fields = meta.fields;
      if (meta.checklistSelection) state.checklistSelection = meta.checklistSelection;
      if (meta.bgStem) state.bgStem = meta.bgStem;
      // 旧バージョン（日付ごとに設定を持っていた形式）からの簡易移行:
      // 代表的な1日付分の設定・BGを、ツール全体の共通設定として引き継ぐ。
      if (!meta.fields && meta.dateSettings) {
        const firstKey = Object.keys(meta.dateSettings)[0];
        if (firstKey) {
          const ds = meta.dateSettings[firstKey];
          if (ds.fields) state.fields = ds.fields;
          if (ds.bgStem) state.bgStem = ds.bgStem;
          if (ds.checklistSelection) state.checklistSelection = ds.checklistSelection;
        }
      }
      if (!meta.selectedDates && meta.currentDate) {
        state.selectedDates = new Set([meta.currentDate]);
      }
    }
    renderAll();
    refreshSidebarTrashCount();
  }

  /* ---------- 2b. ゴミ箱 ----------
     「2. 各ファイルの一覧・ラベル」の各行のゴミ箱アイコンを押すと、そのファイルを読み込んだ
     ファイルの一覧から取り除き、ゴミ箱（IndexedDBの trash ストア）に移す。ゴミ箱の中身の
     閲覧・復元・完全削除は別ページ（trash.html）で行う。 */

  async function refreshSidebarTrashCount() {
    const trashRecords = await dbGetAll(STORE_TRASH);
    const el = document.getElementById('sidebarTrashCount');
    if (el) el.textContent = trashRecords.length > 0 ? `(${trashRecords.length})` : '';
  }

  async function moveStemsToTrash(stems) {
    if (stems.length === 0) { alert('ゴミ箱へ移動するファイルを選択してください。'); return; }
    if (!confirm(`${stems.length}件のファイルをゴミ箱へ移動します。続けますか？（ゴミ箱ページから元に戻せます）`)) return;

    const affectedOthers = new Set();
    for (const stem of stems) {
      if (!state.files[stem]) continue;
      const record = fileEntryToRecord(stem);
      record.trashedAt = Date.now();
      await dbPut(STORE_TRASH, record);
      await dbDeleteKey(STORE_FILES, stem);
      delete state.files[stem];

      if (state.bgStem === stem) state.bgStem = null;
      delete state.checklistSelection.multiFileChecks[stem];
      delete state.checklistSelection.intFileChecks[stem];
      Object.keys(state.files).forEach(other => {
        if (state.files[other].bgOverrideStem === stem) {
          state.files[other].bgOverrideStem = null;
          affectedOthers.add(other);
        }
      });
    }
    affectedOthers.forEach(other => persistFileEntry(other));
    schedulePersistMeta();
    renderAll();
    refreshSidebarTrashCount();
  }

  async function handleFiles(fileList) {
    const arr = Array.from(fileList);
    const touchedStems = new Set();
    for (const file of arr) {
      const ext = extOf(file.name);
      if (ext !== 'csv' && ext !== '887') continue;
      const stem = stemOf(file.name);
      touchedStems.add(stem);
      if (!state.files[stem]) {
        state.files[stem] = {
          name: stem, csvFile: null, condFile: null,
          csvData: null, condData: null,
          label: '', labelType: null, labelValue: null,
          useInPlot: true, dateKey: extractDateKey(stem),
          bgOverrideStem: null,
        };
        registerDateKey(state.files[stem].dateKey);
      }
      if (ext === 'csv') state.files[stem].csvFile = file;
      else state.files[stem].condFile = file;
    }
    for (const stem of touchedStems) {
      const entry = state.files[stem];
      if (entry.condFile && !entry.condData) {
        try {
          const text = await readFileAsText(entry.condFile);
          entry.condData = parseConditions(text);
        } catch (e) {
          console.error('887ファイルの読込に失敗しました', stem, e);
        }
      }
      if (entry.csvFile && !entry.csvData) {
        try {
          await ensureCsvParsed(stem); // 次回自動復元できるよう、読み込み時点で測定データもパースしておく
        } catch (e) {
          console.error('csvファイルの読込に失敗しました', stem, e);
        }
      }
      persistFileEntry(stem);
    }
    schedulePersistMeta();
    renderAll();
  }

  /* ---------- 3. 全体描画 ---------- */

  function cssEscape(s) { return String(s).replace(/[^a-zA-Z0-9_-]/g, '_'); }

  function renderAll() {
    renderSidebarDateList();
    renderFilesTable();
    populateSelects();
  }

  /* ---------- 4. サイドバーの日付一覧（これまで読み込んだ全日付の履歴。複数選択可） ---------- */

  function renderSidebarDateList() {
    const listEl = document.getElementById('sidebarDateList');
    const dateKeys = [...state.knownDates].sort(naturalCompare);
    listEl.innerHTML = '';
    if (dateKeys.length === 0) {
      const p = document.createElement('p');
      p.className = 'sidebar-date-empty';
      p.textContent = 'まだデータがありません';
      listEl.appendChild(p);
      return;
    }
    for (const key of dateKeys) {
      const count = Object.values(state.files).filter(e => e.dateKey === key).length;
      const item = document.createElement('button');
      item.type = 'button';
      const active = state.selectedDates.has(key);
      item.className = 'sidebar-date-item' + (active ? ' active' : '');
      item.setAttribute('aria-pressed', active ? 'true' : 'false');
      item.textContent = `${key} (${count}件)`;
      item.addEventListener('click', () => {
        if (state.selectedDates.has(key)) state.selectedDates.delete(key);
        else state.selectedDates.add(key);
        schedulePersistMeta();
        renderSidebarDateList();
        renderFilesTable();
      });
      listEl.appendChild(item);
    }
  }

  /* ---------- 5. 「2. 各ファイルの一覧・ラベル」（統合表） ---------- */

  function describeLabelType(e) {
    if (e.labelType === 'number') return '数値';
    if (e.labelType === 'string') return '文字列';
    return '(未設定)';
  }

  function updateLabel(stem, value) {
    const e = state.files[stem];
    if (!e) return;
    e.label = value;
    const trimmed = value.trim();
    if (trimmed !== '' && Number.isFinite(Number(trimmed))) {
      e.labelType = 'number';
      e.labelValue = Number(trimmed);
    } else if (trimmed !== '') {
      e.labelType = 'string';
      e.labelValue = trimmed;
    } else {
      e.labelType = null;
      e.labelValue = null;
    }
    const typeEl = document.getElementById(`labelType-${cssEscape(stem)}`);
    if (typeEl) typeEl.textContent = describeLabelType(e);
    populateChecklists();
    persistFileEntry(stem);
  }

  // 個別BGモードON時、ファイルに個別BGの指定があればそれを優先し、なければ共通BGにフォールバックする
  function resolveBgForFile(stem, commonBgStem) {
    const perFileMode = document.getElementById('bgPerFileMode').checked;
    const entry = state.files[stem];
    if (perFileMode && entry && entry.bgOverrideStem && state.files[entry.bgOverrideStem]) {
      return entry.bgOverrideStem;
    }
    return commonBgStem;
  }

  function renderFilesTable() {
    const perFileMode = document.getElementById('bgPerFileMode').checked;
    const allStems = Object.keys(state.files)
      .filter(s => state.selectedDates.has(state.files[s].dateKey))
      .sort((a, b) => {
        const da = state.files[a].dateKey, db = state.files[b].dateKey;
        if (da !== db) return naturalCompare(da, db);
        return naturalCompare(a, b);
      });

    const tbody = document.getElementById('filesTableBody');
    const table = document.getElementById('filesTable');
    const empty = document.getElementById('filesListEmpty');
    document.querySelector('.col-individual-bg').classList.toggle('hidden', !perFileMode);

    if (allStems.length === 0) { table.classList.add('hidden'); empty.classList.remove('hidden'); return; }
    table.classList.remove('hidden'); empty.classList.add('hidden');
    tbody.innerHTML = '';

    // すべての選択済みファイル（個別BGの候補は、表示中の日付に関わらず読み込んだ全ファイルから選べる）
    const allLoadedStems = Object.keys(state.files).sort(naturalCompare);

    let lastDateKey = null;
    for (const stem of allStems) {
      const e = state.files[stem];
      if (e.dateKey !== lastDateKey) {
        lastDateKey = e.dateKey;
        const groupTr = document.createElement('tr');
        groupTr.className = 'date-group-row';
        const groupTd = document.createElement('td');
        groupTd.colSpan = perFileMode ? 12 : 11;
        groupTd.textContent = `日付: ${lastDateKey}`;
        groupTr.appendChild(groupTd);
        tbody.appendChild(groupTr);
      }

      const cond = e.condData || {};
      const tr = document.createElement('tr');

      const tdTrash = document.createElement('td');
      const trashBtn = document.createElement('button');
      trashBtn.type = 'button';
      trashBtn.className = 'trash-row-btn';
      trashBtn.textContent = '🗑';
      trashBtn.title = 'このファイルをゴミ箱へ移動';
      trashBtn.addEventListener('click', () => moveStemsToTrash([stem]));
      tdTrash.appendChild(trashBtn);

      const tdName = document.createElement('td');
      tdName.textContent = stem;
      const tdCond = document.createElement('td');
      tdCond.textContent = e.condFile ? '✓' : '✗';
      const tdCsv = document.createElement('td');
      tdCsv.textContent = e.csvFile ? '✓' : '✗';
      const tdTime = document.createElement('td');
      tdTime.textContent = cond.measurement_time || '';
      const tdSweeps = document.createElement('td');
      tdSweeps.textContent = cond.SWEEPS || '';
      const tdCalfact = document.createElement('td');
      tdCalfact.textContent = ('calfact' in cond) ? formatCalfact(cond.calfact) : '-';

      const tdLabel = document.createElement('td');
      const input = document.createElement('input');
      input.type = 'text';
      input.value = e.label || '';
      input.placeholder = '例: 12.3 または high_power';
      input.addEventListener('input', () => updateLabel(stem, input.value));
      tdLabel.appendChild(input);

      const tdType = document.createElement('td');
      tdType.id = `labelType-${cssEscape(stem)}`;
      tdType.textContent = describeLabelType(e);

      const tdUse = document.createElement('td');
      const useCb = document.createElement('input');
      useCb.type = 'checkbox';
      useCb.className = 'use-in-plot-cb';
      useCb.checked = e.useInPlot !== false;
      useCb.addEventListener('change', () => {
        e.useInPlot = useCb.checked;
        populateChecklists();
        persistFileEntry(stem);
        const useAllCb = document.getElementById('useInPlotSelectAll');
        if (useAllCb) useAllCb.checked = allStems.length > 0 && allStems.every(s => state.files[s].useInPlot !== false);
      });
      tdUse.appendChild(useCb);

      const tdBg = document.createElement('td');
      const radio = document.createElement('input');
      radio.type = 'radio';
      radio.name = 'bgRadio';
      radio.checked = state.bgStem === stem;
      radio.addEventListener('change', () => {
        state.bgStem = stem;
        populateSelects();
        // グラフのBG選択欄を、指定した共通BGファイルに強制的に合わせる
        ['diffBgSelect', 'multiBgSelect', 'intBgSelect'].forEach(id => {
          const sel = document.getElementById(id);
          if (Array.from(sel.options).some(o => o.value === stem)) sel.value = stem;
        });
        populateChecklists();
        schedulePersistMeta();
      });
      tdBg.appendChild(radio);

      const tdIndividualBg = document.createElement('td');
      tdIndividualBg.className = 'individual-bg-cell' + (perFileMode ? '' : ' hidden');
      const bgSelect = document.createElement('select');
      const noneOpt = document.createElement('option');
      noneOpt.value = '';
      noneOpt.textContent = '(共通BGを使用)';
      bgSelect.appendChild(noneOpt);
      for (const other of allLoadedStems) {
        if (other === stem) continue;
        const opt = document.createElement('option');
        opt.value = other;
        opt.textContent = other;
        bgSelect.appendChild(opt);
      }
      bgSelect.value = e.bgOverrideStem && state.files[e.bgOverrideStem] ? e.bgOverrideStem : '';
      bgSelect.addEventListener('change', () => {
        e.bgOverrideStem = bgSelect.value || null;
        persistFileEntry(stem);
      });
      tdIndividualBg.appendChild(bgSelect);

      tr.appendChild(tdTrash);
      tr.appendChild(tdName); tr.appendChild(tdCond); tr.appendChild(tdCsv);
      tr.appendChild(tdTime); tr.appendChild(tdSweeps); tr.appendChild(tdCalfact);
      tr.appendChild(tdLabel); tr.appendChild(tdType); tr.appendChild(tdUse);
      tr.appendChild(tdBg); tr.appendChild(tdIndividualBg);
      tbody.appendChild(tr);
    }
    const useAllCb = document.getElementById('useInPlotSelectAll');
    if (useAllCb) useAllCb.checked = allStems.length > 0 && allStems.every(s => state.files[s].useInPlot !== false);
  }

  /* ---------- 6. 「3. グラフを描画する」向けのファイル選択・チェックリスト ----------
     日付を問わず、読み込んだ全ファイルから選択できる（複数日付のファイルを1つの
     グラフに重ね書きできるようにするため）。 */

  function populateSelects() {
    const stems = Object.keys(state.files).filter(s => state.files[s].csvFile).sort(naturalCompare);
    const selects = ['rawFileSelect', 'massFileSelect', 'diffSignalSelect', 'diffBgSelect', 'multiBgSelect', 'intBgSelect'];
    for (const id of selects) {
      const sel = document.getElementById(id);
      const prevValue = sel.value;
      sel.innerHTML = '';
      for (const stem of stems) {
        const opt = document.createElement('option');
        opt.value = stem;
        opt.textContent = stem;
        sel.appendChild(opt);
      }
      if (stems.includes(prevValue)) {
        sel.value = prevValue;
      } else if ((id === 'diffBgSelect' || id === 'multiBgSelect' || id === 'intBgSelect') &&
                 state.bgStem && stems.includes(state.bgStem)) {
        sel.value = state.bgStem;
      }
    }
    populateChecklists();
  }

  // 数値ラベルが付いているかどうか
  function hasNumericLabel(stem) {
    const e = state.files[stem];
    return !!e && e.labelType === 'number' && Number.isFinite(e.labelValue);
  }

  // 重ね書きの並び順の比較関数。数値ラベルありのファイルはラベル値の昇順で先に並べ、
  // ラベルなしのファイルはその後にファイル名（自然順）で並べる。
  function compareForOverlay(a, b) {
    const na = hasNumericLabel(a), nb = hasNumericLabel(b);
    if (na && nb) return state.files[a].labelValue - state.files[b].labelValue;
    if (na && !nb) return -1;
    if (!na && nb) return 1;
    return naturalCompare(a, b);
  }

  // 重ね書きグラフの各線に表示する名前。数値ラベルがあればラベル値、なければファイル名。
  function overlayTraceName(stem) {
    return hasNumericLabel(stem) ? String(state.files[stem].labelValue) : stem;
  }

  // includeNonNumeric=true の場合は数値ラベルなしのファイルもチェックリストに含める
  // （複数ファイル重ね書き用。ラベルなしはファイル名で識別・並び替えする）。
  function buildChecklist(containerId, stems, excludeStem, includeNonNumeric) {
    const container = document.getElementById(containerId);
    const selection = state.checklistSelection[containerId];
    container.innerHTML = '';
    const sorted = [...stems].sort(compareForOverlay);
    for (const stem of sorted) {
      if (stem === excludeStem) continue;
      const label = document.createElement('label');
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.value = stem;
      cb.checked = selection[stem] !== false; // stem単位で記憶した選択状態（既定はチェック済み）
      cb.addEventListener('change', () => { selection[stem] = cb.checked; schedulePersistMeta(); });
      label.appendChild(cb);
      const text = hasNumericLabel(stem) ? `${stem} (${state.files[stem].labelValue})` : `${stem} (ラベルなし)`;
      label.appendChild(document.createTextNode(text));
      container.appendChild(label);
    }
  }

  function populateChecklists() {
    // 複数ファイル重ね書き: 数値ラベルの有無にかかわらず、プロットに使用する全csvファイルを対象にする
    const plotStems = Object.keys(state.files).filter(s =>
      state.files[s].csvFile && state.files[s].useInPlot !== false);
    // 積算・比プロット: 横軸がラベル値そのものなので、数値ラベルのファイルのみ対象にする
    const numericStems = plotStems.filter(s => hasNumericLabel(s));
    buildChecklist('multiFileChecks', plotStems, document.getElementById('multiBgSelect').value, true);
    buildChecklist('intFileChecks', numericStems, document.getElementById('intBgSelect').value, false);
  }

  /* ---------- 7. 表示設定・軸範囲・質量校正入力値の記憶（ページ更新をまたいで保持する） ---------- */

  const PERSIST_FIELD_IDS = [
    'graphWidth', 'graphHeight', 'graphLineColorEnable', 'graphLineColor', 'graphColorPalette', 'graphLineWidth', 'graphLegendInset',
    'rawShowCond', 'rawShowGrid', 'rawLogY', 'rawXAuto', 'rawXMin', 'rawXMax', 'rawYAuto', 'rawYMin', 'rawYMax',
    'massThreshold', 'massDistance', 'massStep', 'massSmoothWin', 'massRangeAuto', 'massRangeMin', 'massRangeMax',
    'massJacobian', 'massXAuto', 'massXMin', 'massXMax', 'massYAuto', 'massYMin', 'massYMax', 'massLogY', 'massHeLines', 'massCalLines',
    'diffXAxis', 'diffCalibA', 'diffCalibB', 'diffSmoothWin', 'diffShowCond', 'diffLogY',
    'diffXAuto', 'diffXMin', 'diffXMax', 'diffYAuto', 'diffYMin', 'diffYMax',
    'multiXAxis', 'multiCalibA', 'multiCalibB', 'multiSmoothWin', 'multiGrid', 'multiLogY',
    'multiXAuto', 'multiXMin', 'multiXMax', 'multiYAuto', 'multiYMin', 'multiYMax',
    'int1Start', 'int1End', 'int2Start', 'int2End', 'intSmoothWin', 'intXLabel', 'intGrid', 'intErrbar', 'intLine', 'intLogY', 'intLogX',
    'bgPerFileMode',
  ];
  // ファイル選択（rawFileSelect等）・BG選択は、ファイル一覧が変わるたびにoptionを作り直すため
  // ここでの自動復元の対象からは含めない（populateSelects側で個別に復元する）。

  function fieldGet(el) {
    if (el.type === 'checkbox') return el.checked;
    return el.value;
  }
  function fieldSet(el, v) {
    if (el.type === 'checkbox') el.checked = !!v; else el.value = v;
  }

  function bindPersistentFields() {
    for (const id of PERSIST_FIELD_IDS) {
      const el = document.getElementById(id);
      if (!el) continue;
      if (Object.prototype.hasOwnProperty.call(state.fields, id)) {
        fieldSet(el, state.fields[id]);
      } else {
        state.fields[id] = fieldGet(el); // まだ記憶がなければHTMLの既定値を記憶しておく
      }
      const handler = () => {
        state.fields[id] = fieldGet(el);
        schedulePersistMeta();
        if (id === 'bgPerFileMode') renderFilesTable();
      };
      el.addEventListener('input', handler);
      el.addEventListener('change', handler);
    }
  }

  /* ---------- 8. 軸レンジ・グラフ共通ヘルパー ---------- */

  function setXAxisRange(layout, auto, min, max) {
    if (!auto && Number.isFinite(min) && Number.isFinite(max)) layout.xaxis.range = [min, max];
  }
  function setYAxisRange(layout, auto, min, max, log) {
    layout.yaxis.type = log ? 'log' : 'linear';
    if (!auto && Number.isFinite(min) && Number.isFinite(max)) {
      layout.yaxis.range = log
        ? [Math.log10(Math.max(min, 1e-300)), Math.log10(Math.max(max, 1e-300))]
        : [min, max];
    }
  }
  function wireAutoToggle(autoId, minId, maxId) {
    const auto = document.getElementById(autoId);
    const min = document.getElementById(minId);
    const max = document.getElementById(maxId);
    const update = () => { min.disabled = auto.checked; max.disabled = auto.checked; };
    auto.addEventListener('change', update);
    update();
  }

  function syncCalibFields(a, b) {
    if (a == null) return;
    const aStr = a.toExponential(6);
    const bStr = b.toExponential(6);
    for (const id of ['diffCalibA', 'multiCalibA']) {
      const el = document.getElementById(id);
      if (!el.value.trim()) { el.value = aStr; el.dispatchEvent(new Event('change')); }
    }
    for (const id of ['diffCalibB', 'multiCalibB']) {
      const el = document.getElementById(id);
      if (!el.value.trim()) { el.value = bStr; el.dispatchEvent(new Event('change')); }
    }
  }

  // 「3. グラフを描画する」欄で指定された幅・高さ(px)を取得する
  function getGraphSize() {
    const w = parseInt(document.getElementById('graphWidth').value, 10) || 900;
    const h = parseInt(document.getElementById('graphHeight').value, 10) || 450;
    return { w, h };
  }

  // 「⚙ 表示設定」パネルで選べるカラーパレット（1つのグラフに複数の線がある場合の配色）。
  const COLOR_PALETTES = {
    default: null,
    category10: ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    set1: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628', '#f781bf', '#999999'],
    set2: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'],
    pastel: ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#e5d8bd', '#fddaec'],
    blues: ['#08306b', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef'],
  };

  function getGraphStyle() {
    const enableEl = document.getElementById('graphLineColorEnable');
    const colorEl = document.getElementById('graphLineColor');
    const paletteEl = document.getElementById('graphColorPalette');
    const widthEl = document.getElementById('graphLineWidth');
    const legendInsetEl = document.getElementById('graphLegendInset');
    const colorEnabled = enableEl ? enableEl.checked : false;
    return {
      color: colorEnabled && colorEl ? colorEl.value : null,
      palette: paletteEl ? COLOR_PALETTES[paletteEl.value] : null,
      lineWidth: widthEl ? (parseFloat(widthEl.value) || 1) : 1,
      legendInset: legendInsetEl ? legendInsetEl.checked : false,
    };
  }

  function pickColor(palette, index, fallback) {
    if (palette && palette.length) return palette[index % palette.length];
    return fallback;
  }

  function wrapCommaList(text, maxCharsPerLine) {
    if (!text) return text;
    const parts = text.split(', ');
    const lines = [];
    let current = '';
    for (const part of parts) {
      const candidate = current ? `${current}, ${part}` : part;
      if (candidate.length > maxCharsPerLine && current) {
        lines.push(current);
        current = part;
      } else {
        current = candidate;
      }
    }
    if (current) lines.push(current);
    return lines.join('<br>');
  }

  function withMeta(mainTitle, filesText, smoothWindow) {
    const parts = [];
    if (filesText) {
      const { w } = getGraphSize();
      const maxChars = Math.max(30, Math.floor(w / 6.5) - 10);
      parts.push(`File: ${wrapCommaList(filesText, maxChars)}`);
    }
    if (smoothWindow !== undefined && smoothWindow !== null) parts.push(`Smoothing: ${smoothWindow}`);
    if (parts.length === 0) return mainTitle;
    return `${mainTitle}<br><span style="font-size:11px;color:#666666">${parts.join(' &nbsp;|&nbsp; ')}</span>`;
  }

  function sanitizeFilename(s) {
    return String(s).replace(/[^a-zA-Z0-9_\-]+/g, '_');
  }

  function applyLegendLayout(layout) {
    const { legendInset } = getGraphStyle();
    if (legendInset) {
      layout.legend = {
        x: 0.99, xanchor: 'right', y: 0.99, yanchor: 'top', font: { size: 10 },
        bgcolor: 'rgba(255,255,255,0.75)', bordercolor: '#dcdfe4', borderwidth: 1,
      };
    } else {
      if (!layout.margin) layout.margin = {};
      layout.margin.r = Math.max(layout.margin.r || 0, 160);
      layout.legend = { x: 1.02, xanchor: 'left', y: 1, yanchor: 'top', font: { size: 10 } };
    }
  }

  function applyGraphSize(divIds) {
    const { w, h } = getGraphSize();
    divIds.forEach(id => {
      const div = document.getElementById(id);
      if (div) {
        div.style.flex = 'none';
        div.style.width = `${w}px`;
        div.style.height = `${h}px`;
      }
    });
    return { w, h };
  }

  const COPY_ICON = {
    width: 24, height: 24,
    path: 'M19,21H8V7H19M19,5H8A2,2 0 0,0 6,7V21A2,2 0 0,0 8,23H19A2,2 0 0,0 21,21V7A2,2 0 0,0 19,5M16,1H4A2,2 0 0,0 2,3V17H4V3H16V1Z',
  };

  function makeCopyModeBarButton() {
    return {
      name: 'copyImage',
      title: 'グラフをクリップボードにコピー',
      icon: COPY_ICON,
      click: async function (gd) {
        try {
          const url = await Plotly.toImage(gd, {
            format: 'png',
            width: gd._fullLayout.width,
            height: gd._fullLayout.height,
          });
          const blob = await (await fetch(url)).blob();
          await navigator.clipboard.write([new ClipboardItem({ 'image/png': blob })]);
        } catch (e) {
          alert('クリップボードへのコピーに失敗しました。ブラウザでクリップボードへのアクセスが許可されているか確認してください。\n(' + (e && e.message ? e.message : e) + ')');
        }
      },
    };
  }

  function getPlotConfig(filename) {
    return {
      responsive: true,
      modeBarButtonsToAdd: [makeCopyModeBarButton()],
      toImageButtonOptions: { format: 'png', filename: sanitizeFilename(filename || 'graph') },
    };
  }

  /* ---------- 9. タブ切替 ---------- */

  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
      btn.classList.add('active');
      document.getElementById(btn.dataset.tab).classList.add('active');
    });
  });

  /* ---------- 9b. サイドバーのグラフタブへのリンク（クリックでタブ切替＋スクロール） ---------- */
  document.querySelectorAll('.sidebar-graph-tabs a[data-tab-link]').forEach(a => {
    a.addEventListener('click', (ev) => {
      ev.preventDefault();
      const tabId = a.dataset.tabLink;
      const tabBtn = document.querySelector(`.tab-btn[data-tab="${tabId}"]`);
      if (tabBtn) tabBtn.click();
      const target = document.getElementById(tabId);
      if (target) target.scrollIntoView({ behavior: 'smooth', block: 'start' });
    });
  });

  /* ---------- 9c. グラフの「⚙ 表示設定」パネル ---------- */
  const btnGraphSettings = document.getElementById('btnGraphSettings');
  const graphSettingsPanel = document.getElementById('graphSettingsPanel');
  btnGraphSettings.addEventListener('click', () => graphSettingsPanel.classList.toggle('hidden'));
  const graphLineColorEnable = document.getElementById('graphLineColorEnable');
  const graphLineColor = document.getElementById('graphLineColor');
  const syncColorEnabled = () => { graphLineColor.disabled = !graphLineColorEnable.checked; };
  graphLineColorEnable.addEventListener('change', syncColorEnabled);
  syncColorEnabled();

  /* ---------- 10. ファイル入力まわりの配線 ---------- */

  const dropzone = document.getElementById('dropzone');
  const fileInput = document.getElementById('fileInput');
  document.getElementById('btnPickFiles').addEventListener('click', () => fileInput.click());
  fileInput.addEventListener('change', (e) => { handleFiles(e.target.files); fileInput.value = ''; });
  dropzone.addEventListener('dragover', (e) => { e.preventDefault(); dropzone.classList.add('dragover'); });
  dropzone.addEventListener('dragleave', () => dropzone.classList.remove('dragover'));
  dropzone.addEventListener('drop', (e) => {
    e.preventDefault();
    dropzone.classList.remove('dragover');
    handleFiles(e.dataTransfer.files);
  });
  document.getElementById('btnClearFiles').addEventListener('click', () => {
    if (Object.keys(state.files).length > 0 &&
        !confirm('読み込んだファイル・ラベル・表示設定・ゴミ箱の中身はこのブラウザに保存されています。すべてクリアするとゴミ箱の中身も含めて削除され、元に戻せません。続けますか？')) {
      return;
    }
    state.files = {};
    state.knownDates = new Set();
    state.selectedDates = new Set();
    state.bgStem = null;
    state.checklistSelection = { multiFileChecks: {}, intFileChecks: {} };
    dbClear(); // ブラウザに保存していたデータ（ゴミ箱を含む）も合わせて削除する
    document.getElementById('diffCalibA').value = '';
    document.getElementById('diffCalibB').value = '';
    document.getElementById('multiCalibA').value = '';
    document.getElementById('multiCalibB').value = '';
    document.getElementById('massResult').classList.add('hidden');
    document.getElementById('intResultTable').classList.add('hidden');
    ['plotRaw', 'plotMassTof', 'plotMassSpec', 'plotDiffRaw', 'plotDiffDiff',
     'plotMultiRaw', 'plotMultiDiff', 'plotInt1', 'plotInt2', 'plotInt3', 'plotInt4'].forEach(id => {
      const div = document.getElementById(id);
      if (typeof Plotly !== 'undefined' && div.data) Plotly.purge(div);
      div.innerHTML = '';
      delete div.dataset.filename;
    });
    renderAll();
    refreshSidebarTrashCount();
  });

  document.getElementById('useInPlotSelectAll').addEventListener('change', (e) => {
    // ループ中に各行のchangeハンドラがこのヘッダーのcheckedを書き換えてしまうため、
    // 先に「全チェック/全解除」の値を確定させてから各行に反映する（先頭行だけ反映される不具合の対策）。
    const checkAll = e.target.checked;
    const rowCbs = document.querySelectorAll('#filesTableBody tr:not(.date-group-row) .use-in-plot-cb');
    rowCbs.forEach(cb => {
      cb.checked = checkAll;
      cb.dispatchEvent(new Event('change'));
    });
  });

  document.getElementById('multiBgSelect').addEventListener('change', populateChecklists);
  document.getElementById('intBgSelect').addEventListener('change', populateChecklists);

  // 個別BGモード時、Signalファイルを選ぶと、そのファイルに個別BGの指定があれば
  // Backgroundの選択欄に自動でセットする（差分タブでの利便性向上。手動での変更も可能）
  document.getElementById('diffSignalSelect').addEventListener('change', () => {
    const sigStem = document.getElementById('diffSignalSelect').value;
    const entry = state.files[sigStem];
    if (document.getElementById('bgPerFileMode').checked && entry && entry.bgOverrideStem) {
      const bgSel = document.getElementById('diffBgSelect');
      if (Array.from(bgSel.options).some(o => o.value === entry.bgOverrideStem)) {
        bgSel.value = entry.bgOverrideStem;
        bgSel.dispatchEvent(new Event('change'));
      }
    }
  });

  wireAutoToggle('rawXAuto', 'rawXMin', 'rawXMax');
  wireAutoToggle('rawYAuto', 'rawYMin', 'rawYMax');
  wireAutoToggle('massRangeAuto', 'massRangeMin', 'massRangeMax');
  wireAutoToggle('massXAuto', 'massXMin', 'massXMax');
  wireAutoToggle('massYAuto', 'massYMin', 'massYMax');
  wireAutoToggle('diffXAuto', 'diffXMin', 'diffXMax');
  wireAutoToggle('diffYAuto', 'diffYMin', 'diffYMax');
  wireAutoToggle('multiXAuto', 'multiXMin', 'multiXMax');
  wireAutoToggle('multiYAuto', 'multiYMin', 'multiYMax');

  bindPersistentFields();

  /* ============================================================
     11. タブ1: TOFスペクトル（生データ表示）
     ============================================================ */
  document.getElementById('btnPlotRaw').addEventListener('click', async () => {
    const stem = document.getElementById('rawFileSelect').value;
    if (!stem) { alert('ファイルを選択してください'); return; }
    const entry = state.files[stem];
    const csv = await ensureCsvParsed(stem);
    if (!csv) { alert('csvファイルが見つかりません'); return; }
    const cond = entry.condData || {};

    const showCond = document.getElementById('rawShowCond').checked && Object.keys(cond).length > 0;
    const showGrid = document.getElementById('rawShowGrid').checked;
    const logY = document.getElementById('rawLogY').checked;
    const xAuto = document.getElementById('rawXAuto').checked;
    const yAuto = document.getElementById('rawYAuto').checked;
    const xMin = parseFloat(document.getElementById('rawXMin').value);
    const xMax = parseFloat(document.getElementById('rawXMax').value);
    const yMin = parseFloat(document.getElementById('rawYMin').value);
    const yMax = parseFloat(document.getElementById('rawYMax').value);

    applyGraphSize(['plotRaw']);

    const { color: lineColor, lineWidth } = getGraphStyle();
    const trace = { x: csv.time_ns, y: csv.counts, type: 'scattergl', mode: 'lines', line: { width: lineWidth, color: lineColor || 'steelblue' }, name: stem };
    const layout = {
      title: stem,
      xaxis: { title: 'Time (ns)', showgrid: showGrid },
      yaxis: { title: 'Counts', showgrid: showGrid },
      margin: { t: 40 },
    };
    setXAxisRange(layout, xAuto, xMin, xMax);
    setYAxisRange(layout, yAuto, yMin, yMax, logY);
    applyLegendLayout(layout);
    const rawFilename = `TOF_${stem}`;
    Plotly.newPlot('plotRaw', [trace], layout, getPlotConfig(rawFilename));
    document.getElementById('plotRaw').dataset.filename = rawFilename;

    document.getElementById('rawCondText').textContent = showCond ? makeConditionText(cond) : '';
  });

  /* ============================================================
     12. タブ2: 質量校正・質量スペクトル
     ============================================================ */
  document.getElementById('btnPlotMass').addEventListener('click', async () => {
    const stem = document.getElementById('massFileSelect').value;
    if (!stem) { alert('ファイルを選択してください'); return; }
    const entry = state.files[stem];
    const csv = await ensureCsvParsed(stem);
    if (!csv) { alert('csvファイルが見つかりません'); return; }
    const cond = entry.condData || {};
    applyGraphSize(['plotMassTof', 'plotMassSpec']);

    const t_ns = makeTimeAxis(csv, cond);
    const t_s = Float64Array.from(t_ns, v => v * 1e-9);
    const n = t_ns.length;

    const smoothWin = parseInt(document.getElementById('massSmoothWin').value, 10) || 1;
    const countsSmoothed = smooth(csv.counts, smoothWin);

    const threshold = parseFloat(document.getElementById('massThreshold').value) || 0;
    const distance = parseInt(document.getElementById('massDistance').value, 10) || 1;
    let peakIdx = findPeaks(countsSmoothed, threshold, distance);

    const rangeAuto = document.getElementById('massRangeAuto').checked;
    let rangeMsg;
    if (!rangeAuto) {
      const tmin = parseFloat(document.getElementById('massRangeMin').value);
      const tmax = parseFloat(document.getElementById('massRangeMax').value);
      peakIdx = peakIdx.filter(i => t_ns[i] >= tmin && t_ns[i] <= tmax);
      rangeMsg = `検出範囲: ${tmin.toFixed(1)} ～ ${tmax.toFixed(1)} ns`;
    } else {
      rangeMsg = `検出範囲: 全範囲 (${t_ns[0].toFixed(1)} ～ ${t_ns[n - 1].toFixed(1)} ns)`;
    }

    const resultBox = document.getElementById('massResult');
    resultBox.classList.remove('hidden');

    if (peakIdx.length < 2) {
      resultBox.innerHTML = `<span class="warn">エラー: 閾値 ${threshold} counts でピークが ${peakIdx.length} 個しか検出されませんでした。閾値を下げるか、検出範囲・スムージング窓幅を調整してください。</span>`;
      ['plotMassTof', 'plotMassSpec'].forEach(id => {
        const div = document.getElementById(id);
        if (typeof Plotly !== 'undefined' && div.data) Plotly.purge(div);
        div.innerHTML = '';
        delete div.dataset.filename;
      });
      return;
    }

    const massStep = parseFloat(document.getElementById('massStep').value) || 4.0026;
    const peakTimesS = peakIdx.map(i => t_s[i]);
    const calib = calibrate(peakTimesS, massStep);
    const { slope: a, intercept: b, r2, startN } = calib;

    syncCalibFields(a, b);

    resultBox.innerHTML =
      `${rangeMsg}<br>検出されたピーク数: ${peakTimesS.length}<br>` +
      `最適推定: 最初のピークは 質量ステップ×${startN} (質量 ${(massStep * startN).toFixed(2)} u), R²=${r2.toFixed(6)}<br>` +
      `a = ${a.toExponential(4)} s/√u,&nbsp; b = ${b.toExponential(4)} s`;

    const massCalib = new Float64Array(n);
    const jacobian = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const tEff = t_s[i] - b;
      if (tEff > 0) {
        massCalib[i] = (tEff / a) ** 2;
        jacobian[i] = (a * a) / (2 * tEff);
      }
    }
    const useJac = document.getElementById('massJacobian').checked;
    const countsCorrected = useJac ? countsSmoothed.map((v, i) => v * jacobian[i]) : countsSmoothed;
    const intensityLabel = useJac ? 'Intensity (Jacobian Corrected)' : 'Intensity (Smoothed)';
    const plotLabelName = useJac ? 'Jacobian Corrected' : 'Smoothed';

    const { color: lineColor, palette: colorPalette, lineWidth } = getGraphStyle();
    const tofColor = pickColor(colorPalette, 0, lineColor || 'steelblue');
    const peakColor = pickColor(colorPalette, 1, 'red');
    const traceTof = { x: t_ns, y: countsSmoothed, type: 'scattergl', mode: 'lines', line: { width: lineWidth, color: tofColor }, name: 'Smoothed' };
    const tracePeaks = {
      x: peakIdx.map(i => t_ns[i]), y: peakIdx.map(i => countsSmoothed[i]),
      type: 'scatter', mode: 'markers', marker: { color: peakColor, size: 7 }, name: `Peaks (>=${threshold})`,
    };
    const massTofLayout = {
      title: withMeta('TOF Spectrum', stem, smoothWin), xaxis: { title: 'Time (ns)' }, yaxis: { title: 'Counts' }, margin: { t: 40 },
    };
    applyLegendLayout(massTofLayout);
    const massTofFilename = `MassCal_TOF_${stem}`;
    Plotly.newPlot('plotMassTof', [traceTof, tracePeaks], massTofLayout, getPlotConfig(massTofFilename));
    document.getElementById('plotMassTof').dataset.filename = massTofFilename;

    let startIdx = 0;
    while (startIdx < n && t_s[startIdx] <= b) startIdx++;
    const xMass = [], yMass = [];
    for (let i = startIdx; i < n; i++) { xMass.push(massCalib[i]); yMass.push(countsCorrected[i]); }

    const xAuto = document.getElementById('massXAuto').checked;
    const xMinV = parseFloat(document.getElementById('massXMin').value);
    const xMaxV = parseFloat(document.getElementById('massXMax').value);
    const heLines = document.getElementById('massHeLines').checked;
    const calLines = document.getElementById('massCalLines').checked;

    const shapes = [];
    const annotations = [];
    if (heLines) {
      const xmaxForLines = xAuto ? 200 : (Number.isFinite(xMaxV) ? xMaxV : 200);
      for (let m = massStep; m <= xmaxForLines + massStep; m += massStep) {
        shapes.push({ type: 'line', x0: m, x1: m, y0: 0, y1: 1, yref: 'paper', line: { color: 'navy', dash: 'dash', width: 0.5 }, opacity: 0.3 });
      }
    }
    if (calLines) {
      peakTimesS.forEach((tpt, i) => {
        const estMass = massStep * (startN + i);
        shapes.push({ type: 'line', x0: estMass, x1: estMass, y0: 0, y1: 1, yref: 'paper', line: { color: 'green', width: 1 }, opacity: 0.5 });
        annotations.push({ x: estMass, y: 0.85, yref: 'paper', text: `n=${startN + i}`, showarrow: false, font: { color: 'green', size: 9 } });
      });
    }

    const layout2 = {
      title: withMeta(stem, null, smoothWin),
      xaxis: { title: 'Mass / Charge (u)' },
      yaxis: { title: intensityLabel },
      shapes, annotations,
      margin: { t: 40 },
    };
    setXAxisRange(layout2, xAuto, xMinV, xMaxV);
    const yAuto = document.getElementById('massYAuto').checked;
    const yMinV = parseFloat(document.getElementById('massYMin').value);
    const yMaxV = parseFloat(document.getElementById('massYMax').value);
    const logY = document.getElementById('massLogY').checked;
    setYAxisRange(layout2, yAuto, yMinV, yMaxV, logY);
    applyLegendLayout(layout2);

    const massSpecFilename = `MassSpectrum_${stem}`;
    Plotly.newPlot('plotMassSpec', [{ x: xMass, y: yMass, type: 'scattergl', mode: 'lines', line: { width: lineWidth, color: lineColor || 'crimson' }, name: plotLabelName }], layout2, getPlotConfig(massSpecFilename));
    document.getElementById('plotMassSpec').dataset.filename = massSpecFilename;
  });

  /* ============================================================
     13. タブ3: 差分 (Signal - Background)
     ============================================================ */
  document.getElementById('btnPlotDiff').addEventListener('click', async () => {
    const sigStem = document.getElementById('diffSignalSelect').value;
    const bgStem = document.getElementById('diffBgSelect').value;
    if (!sigStem || !bgStem) { alert('SignalとBackgroundファイルを選択してください'); return; }
    if (sigStem === bgStem) { alert('SignalとBackgroundは異なるファイルを選択してください'); return; }

    const sigEntry = state.files[sigStem];
    const bgEntry = state.files[bgStem];
    const sigCsv = await ensureCsvParsed(sigStem);
    const bgCsv = await ensureCsvParsed(bgStem);
    if (!sigCsv || !bgCsv) { alert('csvファイルが見つかりません'); return; }
    applyGraphSize(['plotDiffRaw', 'plotDiffDiff']);
    const sigCond = sigEntry.condData || {};
    const bgCond = bgEntry.condData || {};
    const sigTime = makeTimeAxis(sigCsv, sigCond);
    const bgTime = makeTimeAxis(bgCsv, bgCond);

    const window_ = Math.max(1, parseInt(document.getElementById('diffSmoothWin').value, 10) || 1);
    const sigCounts = smooth(sigCsv.counts, window_);
    const bgCounts = smooth(bgCsv.counts, window_);

    let diffCounts;
    if (arraysEqual(sigTime, bgTime)) {
      diffCounts = sigCounts.map((v, i) => v - bgCounts[i]);
    } else {
      const bgInterp = interpLinear(sigTime, bgTime, bgCounts);
      diffCounts = sigCounts.map((v, i) => v - bgInterp[i]);
    }

    const useMass = document.getElementById('diffXAxis').value === 'mass';
    let xSig, xBgArr, xDiff, xlabel, sigValid, bgValid, diffValid;
    if (useMass) {
      const a = parseFloat(document.getElementById('diffCalibA').value);
      const b = parseFloat(document.getElementById('diffCalibB').value);
      if (!Number.isFinite(a) || !Number.isFinite(b)) {
        alert('質量軸を使う場合は校正パラメータ(a, b)を入力してください。質量校正タブで計算すると自動入力されます。');
        return;
      }
      const sigTS = Float64Array.from(sigTime, v => v * 1e-9);
      const bgTS = Float64Array.from(bgTime, v => v * 1e-9);
      const sigM = timeToMass(sigTS, a, b);
      const bgM = timeToMass(bgTS, a, b);
      const diffM = timeToMass(sigTS, a, b);
      xSig = sigM.mass; sigValid = sigM.mask;
      xBgArr = bgM.mass; bgValid = bgM.mask;
      xDiff = diffM.mass; diffValid = diffM.mask;
      xlabel = 'Mass / Charge (u)';
    } else {
      xSig = sigTime; xBgArr = bgTime; xDiff = sigTime;
      xlabel = 'Time (ns)';
      sigValid = new Uint8Array(sigTime.length).fill(1);
      bgValid = new Uint8Array(bgTime.length).fill(1);
      diffValid = new Uint8Array(sigTime.length).fill(1);
    }

    const sigSeries = filterByMaskPair(xSig, sigCounts, sigValid);
    const bgSeries = filterByMaskPair(xBgArr, bgCounts, bgValid);
    const diffSeries = filterByMaskPair(xDiff, diffCounts, diffValid);

    const xAuto = document.getElementById('diffXAuto').checked;
    const xMin = parseFloat(document.getElementById('diffXMin').value);
    const xMax = parseFloat(document.getElementById('diffXMax').value);
    const yAuto = document.getElementById('diffYAuto').checked;
    const yMin = parseFloat(document.getElementById('diffYMin').value);
    const yMax = parseFloat(document.getElementById('diffYMax').value);
    const logY = document.getElementById('diffLogY').checked;

    const { color: diffLineColor, palette: diffPalette, lineWidth: diffLineWidth } = getGraphStyle();
    const diffSigColor = pickColor(diffPalette, 0, diffLineColor || 'steelblue');
    const diffBgColor = pickColor(diffPalette, 1, 'tomato');

    const layoutRaw = { title: withMeta('Raw Data', `Signal: ${sigStem}, BG: ${bgStem}`, window_), xaxis: { title: xlabel }, yaxis: { title: 'Counts' }, margin: { t: 40 } };
    setXAxisRange(layoutRaw, xAuto, xMin, xMax);
    setYAxisRange(layoutRaw, yAuto, yMin, yMax, logY);
    applyLegendLayout(layoutRaw);
    const diffRawFilename = `Diff_Raw_${sigStem}_vs_${bgStem}`;
    Plotly.newPlot('plotDiffRaw', [
      { x: sigSeries.fx, y: sigSeries.fy, type: 'scattergl', mode: 'lines', name: `Signal: ${sigStem}`, line: { width: diffLineWidth, color: diffSigColor } },
      { x: bgSeries.fx, y: bgSeries.fy, type: 'scattergl', mode: 'lines', name: `BG: ${bgStem}`, line: { width: diffLineWidth, color: diffBgColor }, opacity: 0.7 },
    ], layoutRaw, getPlotConfig(diffRawFilename));
    document.getElementById('plotDiffRaw').dataset.filename = diffRawFilename;

    const layoutDiff = { title: withMeta(`Difference: ${sigStem} − ${bgStem}`, null, window_), xaxis: { title: xlabel }, yaxis: { title: 'Counts (Signal − BG)' }, margin: { t: 40 } };
    setXAxisRange(layoutDiff, xAuto, xMin, xMax);
    setYAxisRange(layoutDiff, yAuto, yMin, yMax, logY);
    applyLegendLayout(layoutDiff);
    const diffDiffFilename = `Diff_${sigStem}_minus_${bgStem}`;
    Plotly.newPlot('plotDiffDiff', [
      { x: diffSeries.fx, y: diffSeries.fy, type: 'scattergl', mode: 'lines', line: { width: diffLineWidth, color: diffLineColor || 'steelblue' }, name: 'diff' },
    ], layoutDiff, getPlotConfig(diffDiffFilename));
    document.getElementById('plotDiffDiff').dataset.filename = diffDiffFilename;

    const showCond = document.getElementById('diffShowCond').checked;
    const condBox = document.getElementById('diffCondText');
    if (showCond) {
      const lines = [];
      lines.push(`Signal: ${sigStem}`);
      lines.push(`  SWEEPS: ${sigCond.SWEEPS || 'N/A'}`);
      if (sigCond.measurement_time) lines.push(`  ${sigCond.measurement_time.slice(0, 45)}`);
      lines.push(`Background: ${bgStem}`);
      lines.push(`  SWEEPS: ${bgCond.SWEEPS || 'N/A'}`);
      if (bgCond.measurement_time) lines.push(`  ${bgCond.measurement_time.slice(0, 45)}`);
      lines.push(`range: ${sigCond.range || 'N/A'}`);
      lines.push(`calfact: ${sigCond.calfact ? formatCalfact(sigCond.calfact) : 'N/A'}`);
      lines.push(`calunit: ${sigCond.calunit || 'N/A'}`);
      lines.push(window_ > 1 ? `smoothing: ${window_}` : 'smoothing: none');
      condBox.textContent = lines.join('\n');
    } else {
      condBox.textContent = '';
    }
  });

  /* ============================================================
     14. タブ4: 複数ファイルの同時プロット (Signal群 - 共通BG / 個別BG)
     複数の日付にまたがるファイルを選んでも、1つのグラフに重ね書きされる。
     ============================================================ */
  document.getElementById('btnPlotMulti').addEventListener('click', async () => {
    const commonBgStem = document.getElementById('multiBgSelect').value;
    if (!commonBgStem) { alert('Backgroundファイルを選択してください'); return; }
    const checked = Array.from(document.querySelectorAll('#multiFileChecks input:checked')).map(c => c.value);
    if (checked.length === 0) { alert('重ね書きする信号ファイルを1つ以上選択してください'); return; }

    applyGraphSize(['plotMultiRaw', 'plotMultiDiff']);
    const window_ = Math.max(1, parseInt(document.getElementById('multiSmoothWin').value, 10) || 1);

    const useMass = document.getElementById('multiXAxis').value === 'mass';
    let calibA, calibB;
    if (useMass) {
      calibA = parseFloat(document.getElementById('multiCalibA').value);
      calibB = parseFloat(document.getElementById('multiCalibB').value);
      if (!Number.isFinite(calibA) || !Number.isFinite(calibB)) {
        alert('質量軸を使う場合は校正パラメータ(a, b)を入力してください');
        return;
      }
    }

    // BGファイルは複数の異なるファイルが使われる可能性があるため、一度読み込んだBGはキャッシュする
    const bgCache = new Map();
    async function getBgSeries(bgStem) {
      if (bgCache.has(bgStem)) return bgCache.get(bgStem);
      const bgEntry = state.files[bgStem];
      const bgCsv = await ensureCsvParsed(bgStem);
      if (!bgCsv) { bgCache.set(bgStem, null); return null; }
      const bgCond = bgEntry.condData || {};
      const bgTime = makeTimeAxis(bgCsv, bgCond);
      const bgCounts = smooth(bgCsv.counts, window_);
      const result = { bgTime, bgCounts };
      bgCache.set(bgStem, result);
      return result;
    }

    const commonBg = await getBgSeries(commonBgStem);
    if (!commonBg) { alert('BGのcsvファイルが見つかりません'); return; }

    const sorted = [...checked].sort(compareForOverlay);

    const rawTraces = [];
    const diffTraces = [];
    let globalXMin = Infinity, globalXMax = -Infinity;
    const { palette: multiPalette, lineWidth: multiLineWidth } = getGraphStyle();
    let fileColorIdx = 0;
    const usedBgStems = new Set();

    for (const stem of sorted) {
      const bgStemForThis = resolveBgForFile(stem, commonBgStem);
      const bgData = await getBgSeries(bgStemForThis);
      if (!bgData) { console.warn(`${bgStemForThis}.csv (BG) が見つからないため ${stem} をスキップします`); continue; }

      const entry = state.files[stem];
      const csv = await ensureCsvParsed(stem);
      if (!csv) { console.warn(`${stem}.csv が見つからないためスキップします`); continue; }
      usedBgStems.add(bgStemForThis);
      const cond = entry.condData || {};
      const sigTime = makeTimeAxis(csv, cond);
      const sigCounts = smooth(csv.counts, window_);
      let diffCounts;
      if (arraysEqual(sigTime, bgData.bgTime)) {
        diffCounts = sigCounts.map((v, i) => v - bgData.bgCounts[i]);
      } else {
        const bgInterp = interpLinear(sigTime, bgData.bgTime, bgData.bgCounts);
        diffCounts = sigCounts.map((v, i) => v - bgInterp[i]);
      }
      const labelText = overlayTraceName(stem);

      let xSig, xDiff, sMask, dMask;
      if (useMass) {
        const r1 = timeToMass(Float64Array.from(sigTime, v => v * 1e-9), calibA, calibB);
        xSig = r1.mass; sMask = r1.mask;
        xDiff = r1.mass; dMask = r1.mask;
      } else {
        xSig = sigTime; xDiff = sigTime;
        sMask = new Uint8Array(sigTime.length).fill(1);
        dMask = sMask;
      }
      const sFiltered = filterByMaskPair(xSig, sigCounts, sMask);
      const dFiltered = filterByMaskPair(xDiff, diffCounts, dMask);
      if (dFiltered.fx.length) {
        globalXMin = Math.min(globalXMin, Math.min(...dFiltered.fx));
        globalXMax = Math.max(globalXMax, Math.max(...dFiltered.fx));
      }
      const fileColor = pickColor(multiPalette, fileColorIdx, undefined);
      fileColorIdx++;
      const rawLine = { width: multiLineWidth };
      const diffLine = { width: multiLineWidth };
      if (fileColor) { rawLine.color = fileColor; diffLine.color = fileColor; }
      rawTraces.push({ x: sFiltered.fx, y: sFiltered.fy, type: 'scattergl', mode: 'lines', name: labelText, line: rawLine });
      diffTraces.push({ x: dFiltered.fx, y: dFiltered.fy, type: 'scattergl', mode: 'lines', name: labelText, line: diffLine });
    }

    // 全ファイルが同じBGを使っている場合のみ、参考として点線のBGトレースを重ねる
    // （個別BGが混在する場合は誤解を招くため表示しない）
    if (usedBgStems.size === 1) {
      const onlyBgStem = [...usedBgStems][0];
      const bgData = await getBgSeries(onlyBgStem);
      if (useMass) {
        const r = timeToMass(Float64Array.from(bgData.bgTime, v => v * 1e-9), calibA, calibB);
        const bgFiltered = filterByMaskPair(r.mass, bgData.bgCounts, r.mask);
        rawTraces.push({ x: bgFiltered.fx, y: bgFiltered.fy, type: 'scattergl', mode: 'lines', name: `BG: ${onlyBgStem}`, line: { width: multiLineWidth, color: 'gray', dash: 'dash' }, opacity: 0.6 });
      } else {
        rawTraces.push({ x: Array.from(bgData.bgTime), y: Array.from(bgData.bgCounts), type: 'scattergl', mode: 'lines', name: `BG: ${onlyBgStem}`, line: { width: multiLineWidth, color: 'gray', dash: 'dash' }, opacity: 0.6 });
      }
    }

    const xlabel = useMass ? 'Mass / Charge (u)' : 'Time (ns)';
    const xAuto = document.getElementById('multiXAuto').checked;
    const xMin = parseFloat(document.getElementById('multiXMin').value);
    const xMax = parseFloat(document.getElementById('multiXMax').value);
    const yAuto = document.getElementById('multiYAuto').checked;
    const yMin = parseFloat(document.getElementById('multiYMin').value);
    const yMax = parseFloat(document.getElementById('multiYMax').value);
    const logY = document.getElementById('multiLogY').checked;
    const showGrid = document.getElementById('multiGrid').checked;

    const bgSummary = usedBgStems.size === 1 ? [...usedBgStems][0] : `${usedBgStems.size}種類（個別BG使用）`;
    const filesText = `${sorted.join(', ')} | BG: ${bgSummary}`;
    const layoutRaw = { title: withMeta('Raw Data', filesText, window_), xaxis: { title: xlabel, showgrid: showGrid }, yaxis: { title: 'Counts', showgrid: showGrid }, margin: { t: 40 } };
    const layoutDiff = { title: withMeta('Difference Plot', filesText, window_), xaxis: { title: xlabel, showgrid: showGrid }, yaxis: { title: 'Counts (Signal − BG)', showgrid: showGrid }, margin: { t: 40 } };
    const xr = xAuto ? [globalXMin, globalXMax] : [xMin, xMax];
    if (Number.isFinite(xr[0]) && Number.isFinite(xr[1])) {
      layoutRaw.xaxis.range = xr; layoutDiff.xaxis.range = xr;
    }
    setYAxisRange(layoutRaw, yAuto, yMin, yMax, logY);
    setYAxisRange(layoutDiff, yAuto, yMin, yMax, logY);
    applyLegendLayout(layoutRaw);
    applyLegendLayout(layoutDiff);

    const bgTag = bgSummary.replace(/[^a-zA-Z0-9_-]+/g, '_');
    const multiRawFilename = `Multi_Raw_BG_${bgTag}`;
    Plotly.newPlot('plotMultiRaw', rawTraces, layoutRaw, getPlotConfig(multiRawFilename));
    document.getElementById('plotMultiRaw').dataset.filename = multiRawFilename;
    const multiDiffFilename = `Multi_Diff_BG_${bgTag}`;
    Plotly.newPlot('plotMultiDiff', diffTraces, layoutDiff, getPlotConfig(multiDiffFilename));
    document.getElementById('plotMultiDiff').dataset.filename = multiDiffFilename;
  });

  /* ============================================================
     15. タブ5: 2範囲積算・比プロット（複数日付をまたいで選択可能）
     ============================================================ */
  document.getElementById('btnPlotIntegration').addEventListener('click', async () => {
    const commonBgStem = document.getElementById('intBgSelect').value;
    if (!commonBgStem) { alert('Backgroundファイルを選択してください'); return; }
    const checked = Array.from(document.querySelectorAll('#intFileChecks input:checked')).map(c => c.value);
    if (checked.length === 0) { alert('積算対象の信号ファイルを1つ以上選択してください'); return; }

    const int1Start = parseFloat(document.getElementById('int1Start').value);
    const int1End = parseFloat(document.getElementById('int1End').value);
    const int2Start = parseFloat(document.getElementById('int2Start').value);
    const int2End = parseFloat(document.getElementById('int2End').value);
    const window_ = Math.max(1, parseInt(document.getElementById('intSmoothWin').value, 10) || 1);
    const xlabel = document.getElementById('intXLabel').value.trim() || 'X';

    applyGraphSize(['plotInt1', 'plotInt2', 'plotInt3', 'plotInt4']);

    const bgCache = new Map();
    async function getBgSeries(bgStem) {
      if (bgCache.has(bgStem)) return bgCache.get(bgStem);
      const bgEntry = state.files[bgStem];
      const bgCsv = await ensureCsvParsed(bgStem);
      if (!bgCsv) { bgCache.set(bgStem, null); return null; }
      const bgCond = bgEntry.condData || {};
      const bgTime = makeTimeAxis(bgCsv, bgCond);
      const result = { bgTime, bgCountsRaw: bgCsv.counts, bgCountsSmooth: smooth(bgCsv.counts, window_) };
      bgCache.set(bgStem, result);
      return result;
    }
    const commonBg = await getBgSeries(commonBgStem);
    if (!commonBg) { alert('BGのcsvファイルが見つかりません'); return; }

    const results = [];
    const usedBgStems = new Set();
    for (const stem of checked) {
      const entry = state.files[stem];
      if (entry.labelType !== 'number') continue;
      const bgStemForThis = resolveBgForFile(stem, commonBgStem);
      const bg = await getBgSeries(bgStemForThis);
      if (!bg) continue;
      const csv = await ensureCsvParsed(stem);
      if (!csv) continue;
      usedBgStems.add(bgStemForThis);
      const cond = entry.condData || {};
      const sigTime = makeTimeAxis(csv, cond);
      const sigCountsRaw = csv.counts;
      const sigCountsSmooth = smooth(sigCountsRaw, window_);

      let diffCounts, bgRawInterp;
      if (arraysEqual(sigTime, bg.bgTime)) {
        diffCounts = sigCountsSmooth.map((v, i) => v - bg.bgCountsSmooth[i]);
        bgRawInterp = bg.bgCountsRaw;
      } else {
        const bgSmoothInterp = interpLinear(sigTime, bg.bgTime, bg.bgCountsSmooth);
        diffCounts = sigCountsSmooth.map((v, i) => v - bgSmoothInterp[i]);
        bgRawInterp = interpLinear(sigTime, bg.bgTime, bg.bgCountsRaw);
      }

      const mask1 = maskRange(sigTime, int1Start, int1End);
      const mask2 = maskRange(sigTime, int2Start, int2End);

      const integral1 = sumMasked(diffCounts, mask1);
      const integral2 = sumMasked(diffCounts, mask2);
      const err1 = Math.sqrt(sumMasked(sigCountsRaw, mask1) + sumMasked(bgRawInterp, mask1));
      const err2 = Math.sqrt(sumMasked(sigCountsRaw, mask2) + sumMasked(bgRawInterp, mask2));

      results.push({ stem, label: entry.labelValue, integral1, integral2, err1, err2, bgStem: bgStemForThis });
    }

    if (results.length === 0) { alert('有効なデータが見つかりませんでした（数値ラベルの信号ファイルを選択してください）'); return; }

    results.sort((a, b) => a.label - b.label);

    const xs = results.map(r => r.label);
    const i1 = results.map(r => r.integral1);
    const i2 = results.map(r => r.integral2);
    const e1 = results.map(r => r.err1);
    const e2 = results.map(r => r.err2);
    const ratio = results.map((r, idx) => (i2[idx] !== 0 ? i1[idx] / i2[idx] : NaN));
    const ratioErr = results.map((r, idx) => {
      if (i1[idx] === 0 || i2[idx] === 0 || !Number.isFinite(ratio[idx])) return NaN;
      return Math.abs(ratio[idx]) * Math.sqrt((e1[idx] / i1[idx]) ** 2 + (e2[idx] / i2[idx]) ** 2);
    });

    const showErr = document.getElementById('intErrbar').checked;
    const logY = document.getElementById('intLogY').checked;
    const logX = document.getElementById('intLogX').checked;
    const useLine = document.getElementById('intLine').checked;
    const showGrid = document.getElementById('intGrid').checked;
    const mode = useLine ? 'lines+markers' : 'markers';

    const { palette: intPalette, lineWidth: intLineWidth } = getGraphStyle();
    function makeTrace(y, err, name, color, symbol, paletteIdx) {
      const resolvedColor = pickColor(intPalette, paletteIdx, color);
      const trace = { x: xs, y, type: 'scatter', mode, name, marker: { color: resolvedColor, symbol }, line: { color: resolvedColor, width: intLineWidth } };
      if (showErr) trace.error_y = { type: 'data', array: err, visible: true };
      return trace;
    }
    const bgSummary = usedBgStems.size === 1 ? [...usedBgStems][0] : `${usedBgStems.size}種類（個別BG使用）`;
    const filesText = `${results.map(r => r.stem).join(', ')} | BG: ${bgSummary}`;
    const baseLayout = (title, ytitle) => {
      const layout = {
        title: withMeta(title, filesText, window_),
        xaxis: { title: xlabel, type: logX ? 'log' : 'linear', showgrid: showGrid },
        yaxis: { title: ytitle, type: logY ? 'log' : 'linear', showgrid: showGrid },
        margin: { t: 40 },
      };
      applyLegendLayout(layout);
      return layout;
    };

    const bgTag = bgSummary.replace(/[^a-zA-Z0-9_-]+/g, '_');
    const int1Filename = `Integration_Range1_BG_${bgTag}`;
    Plotly.newPlot('plotInt1', [makeTrace(i1, e1, 'Range 1', 'blue', 'circle', 0)], baseLayout(`Range 1 (${int1Start} - ${int1End} ns)`, 'Integrated Counts'), getPlotConfig(int1Filename));
    document.getElementById('plotInt1').dataset.filename = int1Filename;
    const int2Filename = `Integration_Range2_BG_${bgTag}`;
    Plotly.newPlot('plotInt2', [makeTrace(i2, e2, 'Range 2', 'orange', 'square', 0)], baseLayout(`Range 2 (${int2Start} - ${int2End} ns)`, 'Integrated Counts'), getPlotConfig(int2Filename));
    document.getElementById('plotInt2').dataset.filename = int2Filename;
    const int3Filename = `Integration_Simultaneous_BG_${bgTag}`;
    Plotly.newPlot('plotInt3', [makeTrace(i1, e1, 'Range 1', 'blue', 'circle', 0), makeTrace(i2, e2, 'Range 2', 'orange', 'square', 1)], baseLayout('Simultaneous Plot', 'Integrated Counts'), getPlotConfig(int3Filename));
    document.getElementById('plotInt3').dataset.filename = int3Filename;
    const int4Filename = `Integration_Ratio_BG_${bgTag}`;
    Plotly.newPlot('plotInt4', [makeTrace(ratio, ratioErr, 'Range 1 / Range 2', 'green', 'triangle-up', 0)], baseLayout('Ratio (Range 1 / Range 2)', 'Range 1 / Range 2'), getPlotConfig(int4Filename));
    document.getElementById('plotInt4').dataset.filename = int4Filename;

    const tbody = document.getElementById('intResultTableBody');
    tbody.innerHTML = '';
    results.forEach((r, idx) => {
      const tr = document.createElement('tr');
      const rStr = Number.isFinite(ratio[idx]) ? `${ratio[idx].toFixed(3)} ± ${ratioErr[idx].toFixed(3)}` : 'N/A';
      tr.innerHTML = `<td>${r.stem}</td><td>${r.label}</td><td>${r.integral1.toExponential(3)} ± ${r.err1.toExponential(3)}</td><td>${r.integral2.toExponential(3)} ± ${r.err2.toExponential(3)}</td><td>${rStr}</td>`;
      tbody.appendChild(tr);
    });
    document.getElementById('intResultTable').classList.remove('hidden');
  });

  /* ---------- 初期描画（IndexedDBに保存済みのファイル・設定があれば自動復元する） ---------- */
  renderAll();
  restoreFromDB();
}
