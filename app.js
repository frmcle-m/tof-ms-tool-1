/* ============================================================
   TOF-MS データ解析ツール — 実装コード
   (このファイルの中身は index.html の画面には表示されません)
   ============================================================ */

/* ============================================================
   0. グローバル状態
   ============================================================ */
const state = {
  files: {},          // stem -> { name, csvFile, condFile, csvData, condData, label, labelType, labelValue, dateKey, useInPlot }
  bgStem: null,       // ラベル表で「BGとして使用」に指定されたファイル
  lastCalib: { a: null, b: null }, // 質量校正タブで計算された a, b
  currentDate: null,  // 現在選択中の日付タブ（ファイル名先頭の日付キー）
  checklistSelection: { multiFileChecks: {}, intFileChecks: {} }, // stem -> true/false（チェック状態をstem単位で記憶し、日付切替やDOM再描画で消えないようにする）
};

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
function makeConditionText(cond) {
  const lines = [];
  for (const key in KEY_LABELS) {
    if (key in cond) lines.push(`${KEY_LABELS[key]}: ${cond[key]}`);
  }
  return lines.join('\n');
}

/* ============================================================
   以下はブラウザ環境（DOM）でのみ実行する
   ============================================================ */
if (typeof document !== 'undefined') {

  /* ---------- 2. ファイル読込ヘルパー ---------- */

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

  /* ---------- 2b. ローカル永続化 (IndexedDB) ----------
     読み込んだファイル（測定データ含む）とラベル・設定をブラウザのIndexedDBに保存し、
     次回このHTMLを開いたときに自動で復元する。データは常にこのブラウザ内にだけ保存され、
     外部サーバーには一切送信されない。 */
  const DB_NAME = 'tof_tool_db';
  const DB_VERSION = 1;
  const STORE_FILES = 'files';
  const STORE_META = 'meta';

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

  async function dbClear() {
    const db = await dbPromise;
    if (!db) return;
    return new Promise((resolve) => {
      const tx = db.transaction([STORE_FILES, STORE_META], 'readwrite');
      tx.objectStore(STORE_FILES).clear();
      tx.objectStore(STORE_META).clear();
      tx.oncomplete = () => resolve();
      tx.onerror = () => resolve();
    });
  }

  // 1ファイル分（測定データ・条件・ラベル等）をIndexedDBに保存する
  function persistFileEntry(stem) {
    const e = state.files[stem];
    if (!e) return;
    dbPut(STORE_FILES, {
      stem,
      csvData: e.csvData ? { time_ns: Array.from(e.csvData.time_ns), counts: Array.from(e.csvData.counts) } : null,
      condData: e.condData || {},
      label: e.label, labelType: e.labelType, labelValue: e.labelValue,
      useInPlot: e.useInPlot, dateKey: e.dateKey,
      hasCsv: !!e.csvFile, hasCond: !!e.condFile,
    });
  }

  // BG指定・現在の日付タブなど、ファイル単位ではない設定をIndexedDBに保存する
  function persistMeta() {
    dbPut(STORE_META, {
      key: 'settings',
      bgStem: state.bgStem,
      currentDate: state.currentDate,
      checklistSelection: state.checklistSelection,
    });
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
      };
    }
    const metaRecords = await dbGetAll(STORE_META);
    const meta = metaRecords.find(m => m.key === 'settings');
    if (meta) {
      state.bgStem = meta.bgStem || null;
      state.currentDate = meta.currentDate || null;
      if (meta.checklistSelection) state.checklistSelection = meta.checklistSelection;
    }
    renderAll();
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
        };
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
    renderAll();
  }

  /* ---------- 3. 画面描画（一覧・ラベル） ---------- */

  function cssEscape(s) { return s.replace(/[^a-zA-Z0-9_-]/g, '_'); }

  function renderAll() {
    renderDateTabs();
    renderFileSummary();
    renderFileTable();
    renderLabelTable();
    populateSelects();
  }

  // 読み込み済みファイルを日付キーごとにグループ化し、タブとして表示する。
  // state.currentDate がタブ選択と連動し、他の描画関数はこの値でファイルを絞り込む。
  function renderDateTabs() {
    const wrap = document.getElementById('section-date');
    const container = document.getElementById('dateTabs');
    const dateKeys = [...new Set(Object.values(state.files).map(e => e.dateKey))].sort(naturalCompare);
    if (dateKeys.length === 0) {
      wrap.classList.add('hidden');
      state.currentDate = null;
      container.innerHTML = '';
      return;
    }
    wrap.classList.remove('hidden');
    if (!state.currentDate || !dateKeys.includes(state.currentDate)) {
      state.currentDate = dateKeys[0];
    }
    container.innerHTML = '';
    for (const key of dateKeys) {
      const count = Object.values(state.files).filter(e => e.dateKey === key).length;
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'date-tab-btn' + (key === state.currentDate ? ' active' : '');
      btn.textContent = `${key} (${count}件)`;
      btn.addEventListener('click', () => {
        state.currentDate = key;
        persistMeta();
        renderAll();
      });
      container.appendChild(btn);
    }
  }

  function renderFileSummary() {
    const allStems = Object.keys(state.files);
    const el = document.getElementById('fileSummary');
    if (allStems.length === 0) { el.innerHTML = ''; return; }
    const stems = allStems.filter(s => state.files[s].dateKey === state.currentDate);
    let paired = 0, csvOnly = 0, condOnly = 0;
    for (const s of stems) {
      const e = state.files[s];
      if (e.csvFile && e.condFile) paired++;
      else if (e.csvFile) csvOnly++;
      else condOnly++;
    }
    let html = `<span class="ok">${stems.length} 件のファイル名を認識（csv/887ペア: ${paired} 組）— 日付: ${state.currentDate}</span>`;
    if (csvOnly > 0) html += `<br><span class="warn">.887 が見つからないファイル: ${csvOnly} 件</span>`;
    if (condOnly > 0) html += `<br><span class="warn">.csv が見つからないファイル: ${condOnly} 件</span>`;
    if (allStems.length !== stems.length) {
      html += `<br><span class="hint-inline">（読み込み済みの全ファイル: ${allStems.length} 件。他の日付は上の「日付を選択」タブから見られます）</span>`;
    }
    el.innerHTML = html;
  }

  function renderFileTable() {
    const stems = Object.keys(state.files)
      .filter(s => state.files[s].dateKey === state.currentDate)
      .sort(naturalCompare);
    const tbody = document.getElementById('fileTableBody');
    const table = document.getElementById('fileTable');
    const empty = document.getElementById('fileListEmpty');
    if (stems.length === 0) { table.classList.add('hidden'); empty.classList.remove('hidden'); return; }
    table.classList.remove('hidden'); empty.classList.add('hidden');
    tbody.innerHTML = '';
    for (const stem of stems) {
      const e = state.files[stem];
      const cond = e.condData || {};
      const tr = document.createElement('tr');
      tr.innerHTML = `
        <td>${stem}</td>
        <td>${e.csvFile ? '✓' : '✗'}</td>
        <td>${e.condFile ? '✓' : '✗'}</td>
        <td>${cond.measurement_time || ''}</td>
        <td>${cond.SWEEPS || ''}</td>`;
      tbody.appendChild(tr);
    }
  }

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

  function renderLabelTable() {
    const stems = Object.keys(state.files)
      .filter(s => state.files[s].dateKey === state.currentDate)
      .sort(naturalCompare);
    const tbody = document.getElementById('labelTableBody');
    const table = document.getElementById('labelTable');
    const empty = document.getElementById('labelListEmpty');
    if (stems.length === 0) { table.classList.add('hidden'); empty.classList.remove('hidden'); return; }
    table.classList.remove('hidden'); empty.classList.add('hidden');
    tbody.innerHTML = '';
    for (const stem of stems) {
      const e = state.files[stem];
      const tr = document.createElement('tr');

      const tdName = document.createElement('td');
      tdName.textContent = stem;

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

      const tdCalfact = document.createElement('td');
      const cond = e.condData || {};
      tdCalfact.textContent = ('calfact' in cond) ? cond.calfact : '-';

      const tdUse = document.createElement('td');
      const useCb = document.createElement('input');
      useCb.type = 'checkbox';
      useCb.checked = e.useInPlot !== false;
      useCb.addEventListener('change', () => {
        e.useInPlot = useCb.checked;
        populateChecklists();
        persistFileEntry(stem);
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
        // グラフのBG選択欄を、ラベルタブで新しく指定したBGファイルに強制的に合わせる
        // （populateSelectsは既存の選択を優先して保持するため、ここで明示的に上書きする）
        ['diffBgSelect', 'multiBgSelect', 'intBgSelect'].forEach(id => {
          const sel = document.getElementById(id);
          if (Array.from(sel.options).some(o => o.value === stem)) sel.value = stem;
        });
        populateChecklists();
        persistMeta();
      });
      tdBg.appendChild(radio);

      tr.appendChild(tdName); tr.appendChild(tdLabel); tr.appendChild(tdType);
      tr.appendChild(tdCalfact); tr.appendChild(tdUse); tr.appendChild(tdBg);
      tbody.appendChild(tr);
    }
  }

  function populateSelects() {
    const stems = Object.keys(state.files)
      .filter(s => state.files[s].csvFile && state.files[s].dateKey === state.currentDate)
      .sort(naturalCompare);
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

  function buildChecklist(containerId, stems, excludeStem) {
    const container = document.getElementById(containerId);
    const selection = state.checklistSelection[containerId];
    container.innerHTML = '';
    const sorted = [...stems].sort((a, b) => state.files[a].labelValue - state.files[b].labelValue);
    // naturalCompareでファイル名順に並べたい場合は上の行を
    // const sorted = [...stems].sort(naturalCompare); に置き換えてください（既定はラベル値の昇順）
    for (const stem of sorted) {
      if (stem === excludeStem) continue;
      const label = document.createElement('label');
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.value = stem;
      cb.checked = selection[stem] !== false; // stem単位で記憶した選択状態（既定はチェック済み）
      cb.addEventListener('change', () => { selection[stem] = cb.checked; });
      label.appendChild(cb);
      label.appendChild(document.createTextNode(`${stem} (${state.files[stem].labelValue})`));
      container.appendChild(label);
    }
  }

  function populateChecklists() {
    const numericStems = Object.keys(state.files).filter(s =>
      state.files[s].csvFile && state.files[s].labelType === 'number' &&
      state.files[s].useInPlot !== false && state.files[s].dateKey === state.currentDate);
    buildChecklist('multiFileChecks', numericStems, document.getElementById('multiBgSelect').value);
    buildChecklist('intFileChecks', numericStems, document.getElementById('intBgSelect').value);
  }

  /* ---------- 4. 軸レンジ・グラフ共通ヘルパー ---------- */

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

  function syncCalibFields() {
    if (state.lastCalib.a == null) return;
    const aStr = state.lastCalib.a.toExponential(6);
    const bStr = state.lastCalib.b.toExponential(6);
    for (const id of ['diffCalibA', 'multiCalibA']) {
      const el = document.getElementById(id);
      if (!el.value.trim()) el.value = aStr;
    }
    for (const id of ['diffCalibB', 'multiCalibB']) {
      const el = document.getElementById(id);
      if (!el.value.trim()) el.value = bStr;
    }
  }

  // グラフタイトルに「使用ファイル」「スムージング窓幅」を追記する（PNG保存時にも画像に焼き込まれる）
  // 「4. グラフを描画する」欄で指定された幅・高さ(px)を取得する
  function getGraphSize() {
    const w = parseInt(document.getElementById('graphWidth').value, 10) || 900;
    const h = parseInt(document.getElementById('graphHeight').value, 10) || 450;
    return { w, h };
  }

  // 「⚙ 表示設定」パネルで選べるカラーパレット（1つのグラフに複数の線がある場合の配色）。
  // 「既定」(default) を選んだ場合は null を返し、これまで通りの固定色（各グラフごとの
  // steelblue/tomato/redなど）をそのまま使う。他のパレットを選んだ場合のみ、配列の色を
  // 順番に割り当てて上書きする。
  const COLOR_PALETTES = {
    default: null,
    category10: ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    set1: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628', '#f781bf', '#999999'],
    set2: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'],
    pastel: ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#e5d8bd', '#fddaec'],
    blues: ['#08306b', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef'],
  };

  // 「⚙ 表示設定」パネルの現在の値（データ線の色・カラーパレット・線の太さ・凡例のグラフ内表示）を取得する。
  // 「データ線の色を指定する」がオフの場合、colorはnullを返し、各グラフはこれまで通りの
  // 既定色（steelblue/crimsonなど）をそのまま使う。各入力欄が存在しない場合（読み込み順の
  // 都合など）も、これまでの見た目と同じ既定値を返す。
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

  // 1つのグラフに複数の線がある場合の色を決める。パレットが選択されていればその配列から
  // index番目の色を順に割り当て、「既定」パレットの場合はこれまで通り固定色(fallback)を使う。
  function pickColor(palette, index, fallback) {
    if (palette && palette.length) return palette[index % palette.length];
    return fallback;
  }

  // カンマ区切りの長いリスト（ファイル名など）を、指定した1行あたりの文字数を超えないよう
  // カンマの区切りごとに改行(<br>)を挿入する。ファイル名の途中では改行しない。
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

  // グラフタイトルに「使用ファイル」「スムージング窓幅」を追記する（PNG保存時にも画像に焼き込まれる）。
  // ファイル名の合計が現在のグラフ幅より長くなる場合は自動的に改行する。
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

  // ファイル名に使えない文字を "_" に置き換える
  function sanitizeFilename(s) {
    return String(s).replace(/[^a-zA-Z0-9_\-]+/g, '_');
  }

  // 凡例(レジェンド)の配置。「⚙ 表示設定」で「凡例をグラフ内に表示(inset)」がオフの場合は
  // 従来通りグラフの外・右側に配置し、ファイル数が多くても凡例がグラフ内に収まりきらず
  // 切れてしまうのを防ぐ。オンの場合はグラフ内側（右上）に重ねて表示する。
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

  // 「4. グラフを描画する」欄で指定された幅・高さを、指定した各プロットのコンテナに反映する
  function applyGraphSize(divIds) {
    const { w, h } = getGraphSize();
    divIds.forEach(id => {
      const div = document.getElementById(id);
      if (div) {
        div.style.flex = 'none'; // 幅・高さを固定するため、plot-areaのflex伸縮指定を無効化する
        div.style.width = `${w}px`;
        div.style.height = `${h}px`;
      }
    });
    return { w, h };
  }

  // 各グラフ右上のツールバー（モードバー）に表示する「コピー」ボタンのアイコン。
  // MDI (Material Design Icons) の content-copy と同じ形状（2枚の四角が重なった、
  // 一般的な「コピー」を表すアイコン）。
  const COPY_ICON = {
    width: 24, height: 24,
    path: 'M19,21H8V7H19M19,5H8A2,2 0 0,0 6,7V21A2,2 0 0,0 8,23H19A2,2 0 0,0 21,21V7A2,2 0 0,0 19,5M16,1H4A2,2 0 0,0 2,3V17H4V3H16V1Z',
  };

  // グラフ画像をPNGとしてクリップボードにコピーする、モードバー用のカスタムボタン。
  // 標準のカメラ（ダウンロード）アイコンの隣に表示される。コピーしたPNGはそのまま
  // Word・PowerPoint・Slackなどに貼り付けて使える。
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

  // Plotly.newPlot に渡す共通config。「コピー」ボタンの追加と、カメラ（ダウンロード）
  // アイコンでPNG保存する際の既定ファイル名をまとめて設定する。
  function getPlotConfig(filename) {
    return {
      responsive: true,
      modeBarButtonsToAdd: [makeCopyModeBarButton()],
      toImageButtonOptions: { format: 'png', filename: sanitizeFilename(filename || 'graph') },
    };
  }

  /* ---------- 5. タブ切替 ---------- */

  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
      btn.classList.add('active');
      document.getElementById(btn.dataset.tab).classList.add('active');
    });
  });

  /* ---------- 5b. グラフの「⚙ 表示設定」パネル ----------
     グラフサイズ入力の隣にある歯車ボタンで、色・線の太さ・凡例配置などの
     パネルの表示/非表示を切り替える。設定値は各グラフを「表示」ボタンで
     描画する際に読み込まれる（グラフサイズと同じ扱い）。 */
  const btnGraphSettings = document.getElementById('btnGraphSettings');
  const graphSettingsPanel = document.getElementById('graphSettingsPanel');
  if (btnGraphSettings && graphSettingsPanel) {
    btnGraphSettings.addEventListener('click', () => {
      graphSettingsPanel.classList.toggle('hidden');
    });
  }
  const graphLineColorEnable = document.getElementById('graphLineColorEnable');
  const graphLineColor = document.getElementById('graphLineColor');
  if (graphLineColorEnable && graphLineColor) {
    const syncColorEnabled = () => { graphLineColor.disabled = !graphLineColorEnable.checked; };
    graphLineColorEnable.addEventListener('change', syncColorEnabled);
    syncColorEnabled();
  }

  /* ---------- 6. ファイル入力まわりの配線 ---------- */

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
        !confirm('読み込んだファイル・ラベルはこのブラウザに保存されています。すべてクリアすると保存内容も削除され、元に戻せません。続けますか？')) {
      return;
    }
    state.files = {};
    state.bgStem = null;
    state.currentDate = null;
    state.lastCalib = { a: null, b: null };
    state.checklistSelection = { multiFileChecks: {}, intFileChecks: {} };
    dbClear(); // ブラウザに保存していたデータも合わせて削除する
    document.getElementById('diffCalibA').value = '';
    document.getElementById('diffCalibB').value = '';
    document.getElementById('multiCalibA').value = '';
    document.getElementById('multiCalibB').value = '';
    document.getElementById('massResult').classList.add('hidden');
    document.getElementById('intResultTable').classList.add('hidden');
    ['plotRaw', 'plotMassTof', 'plotMassSpec', 'plotDiffRaw', 'plotDiffDiff',
     'plotMultiRaw', 'plotMultiDiff', 'plotInt1', 'plotInt2', 'plotInt3', 'plotInt4'].forEach(id => {
      const div = document.getElementById(id);
      if (typeof Plotly !== 'undefined' && div.data) Plotly.purge(div); // .data/.layoutなどPlotlyの内部状態も完全に破棄する
      div.innerHTML = '';
      delete div.dataset.filename;
    });
    renderAll();
  });

  document.getElementById('multiBgSelect').addEventListener('change', populateChecklists);
  document.getElementById('intBgSelect').addEventListener('change', populateChecklists);

  wireAutoToggle('rawXAuto', 'rawXMin', 'rawXMax');
  wireAutoToggle('rawYAuto', 'rawYMin', 'rawYMax');
  wireAutoToggle('massRangeAuto', 'massRangeMin', 'massRangeMax');
  wireAutoToggle('massXAuto', 'massXMin', 'massXMax');
  wireAutoToggle('massYAuto', 'massYMin', 'massYMax');
  wireAutoToggle('diffXAuto', 'diffXMin', 'diffXMax');
  wireAutoToggle('diffYAuto', 'diffYMin', 'diffYMax');
  wireAutoToggle('multiXAuto', 'multiXMin', 'multiXMax');
  wireAutoToggle('multiYAuto', 'multiYMin', 'multiYMax');

  /* ============================================================
     7. タブ1: TOFスペクトル（生データ表示）
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
    const rawFilename = `${state.currentDate}_TOF_${stem}`;
    Plotly.newPlot('plotRaw', [trace], layout, getPlotConfig(rawFilename));
    document.getElementById('plotRaw').dataset.filename = rawFilename;

    document.getElementById('rawCondText').textContent = showCond ? makeConditionText(cond) : '';
  });

  /* ============================================================
     8. タブ2: 質量校正・質量スペクトル
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

    state.lastCalib = { a, b };
    syncCalibFields();

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

    // 上段: TOFスペクトル + 検出ピーク（この2系列は「⚙ 表示設定」のパレット選択時、パレットの色を順に使う）
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
    const massTofFilename = `${state.currentDate}_MassCal_TOF_${stem}`;
    Plotly.newPlot('plotMassTof', [traceTof, tracePeaks], massTofLayout, getPlotConfig(massTofFilename));
    document.getElementById('plotMassTof').dataset.filename = massTofFilename;

    // 下段: 質量スペクトル
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

    const massSpecFilename = `${state.currentDate}_MassSpectrum_${stem}`;
    Plotly.newPlot('plotMassSpec', [{ x: xMass, y: yMass, type: 'scattergl', mode: 'lines', line: { width: lineWidth, color: lineColor || 'crimson' }, name: plotLabelName }], layout2, getPlotConfig(massSpecFilename));
    document.getElementById('plotMassSpec').dataset.filename = massSpecFilename;
  });

  /* ============================================================
     9. タブ3: 差分 (Signal - Background)
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
    const diffRawFilename = `${state.currentDate}_Diff_Raw_${sigStem}_vs_${bgStem}`;
    Plotly.newPlot('plotDiffRaw', [
      { x: sigSeries.fx, y: sigSeries.fy, type: 'scattergl', mode: 'lines', name: `Signal: ${sigStem}`, line: { width: diffLineWidth, color: diffSigColor } },
      { x: bgSeries.fx, y: bgSeries.fy, type: 'scattergl', mode: 'lines', name: `BG: ${bgStem}`, line: { width: diffLineWidth, color: diffBgColor }, opacity: 0.7 },
    ], layoutRaw, getPlotConfig(diffRawFilename));
    document.getElementById('plotDiffRaw').dataset.filename = diffRawFilename;

    const layoutDiff = { title: withMeta(`Difference: ${sigStem} − ${bgStem}`, null, window_), xaxis: { title: xlabel }, yaxis: { title: 'Counts (Signal − BG)' }, margin: { t: 40 } };
    setXAxisRange(layoutDiff, xAuto, xMin, xMax);
    setYAxisRange(layoutDiff, yAuto, yMin, yMax, logY);
    applyLegendLayout(layoutDiff);
    const diffDiffFilename = `${state.currentDate}_Diff_${sigStem}_minus_${bgStem}`;
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
      lines.push(`calfact: ${sigCond.calfact || 'N/A'}`);
      lines.push(`calunit: ${sigCond.calunit || 'N/A'}`);
      lines.push(window_ > 1 ? `smoothing: ${window_}` : 'smoothing: none');
      condBox.textContent = lines.join('\n');
    } else {
      condBox.textContent = '';
    }
  });

  /* ============================================================
     10. タブ4: 複数ファイルの同時プロット (Signal群 - 共通BG)
     ============================================================ */
  document.getElementById('btnPlotMulti').addEventListener('click', async () => {
    const bgStem = document.getElementById('multiBgSelect').value;
    if (!bgStem) { alert('Backgroundファイルを選択してください'); return; }
    const checked = Array.from(document.querySelectorAll('#multiFileChecks input:checked')).map(c => c.value);
    if (checked.length === 0) { alert('重ね書きする信号ファイルを1つ以上選択してください'); return; }

    const bgEntry = state.files[bgStem];
    const bgCsv = await ensureCsvParsed(bgStem);
    if (!bgCsv) { alert('BGのcsvファイルが見つかりません'); return; }
    applyGraphSize(['plotMultiRaw', 'plotMultiDiff']);
    const bgCond = bgEntry.condData || {};
    const bgTime = makeTimeAxis(bgCsv, bgCond);

    const window_ = Math.max(1, parseInt(document.getElementById('multiSmoothWin').value, 10) || 1);
    const bgCounts = smooth(bgCsv.counts, window_);

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

    let bgMassX, bgMassMask;
    if (useMass) {
      const r = timeToMass(Float64Array.from(bgTime, v => v * 1e-9), calibA, calibB);
      bgMassX = r.mass; bgMassMask = r.mask;
    }

    const sorted = [...checked].sort((s1, s2) => state.files[s1].labelValue - state.files[s2].labelValue);

    const rawTraces = [];
    const diffTraces = [];
    let globalXMin = Infinity, globalXMax = -Infinity;
    const { palette: multiPalette, lineWidth: multiLineWidth } = getGraphStyle();
    let fileColorIdx = 0;

    for (const stem of sorted) {
      const entry = state.files[stem];
      const csv = await ensureCsvParsed(stem);
      if (!csv) { console.warn(`${stem}.csv が見つからないためスキップします`); continue; }
      const cond = entry.condData || {};
      const sigTime = makeTimeAxis(csv, cond);
      const sigCounts = smooth(csv.counts, window_);
      let diffCounts;
      if (arraysEqual(sigTime, bgTime)) {
        diffCounts = sigCounts.map((v, i) => v - bgCounts[i]);
      } else {
        const bgInterp = interpLinear(sigTime, bgTime, bgCounts);
        diffCounts = sigCounts.map((v, i) => v - bgInterp[i]);
      }
      const labelText = String(entry.labelValue);

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
      // ファイルごとに異なる色を割り当てる。パレットが「既定」の場合はcolorを指定せず、
      // これまで通りPlotlyの自動配色に任せる。
      const fileColor = pickColor(multiPalette, fileColorIdx, undefined);
      fileColorIdx++;
      const rawLine = { width: multiLineWidth };
      const diffLine = { width: multiLineWidth };
      if (fileColor) { rawLine.color = fileColor; diffLine.color = fileColor; }
      rawTraces.push({ x: sFiltered.fx, y: sFiltered.fy, type: 'scattergl', mode: 'lines', name: labelText, line: rawLine });
      diffTraces.push({ x: dFiltered.fx, y: dFiltered.fy, type: 'scattergl', mode: 'lines', name: labelText, line: diffLine });
    }

    if (useMass) {
      const bgFiltered = filterByMaskPair(bgMassX, bgCounts, bgMassMask);
      rawTraces.push({ x: bgFiltered.fx, y: bgFiltered.fy, type: 'scattergl', mode: 'lines', name: `BG: ${bgStem}`, line: { width: multiLineWidth, color: 'gray', dash: 'dash' }, opacity: 0.6 });
    } else {
      rawTraces.push({ x: Array.from(bgTime), y: Array.from(bgCounts), type: 'scattergl', mode: 'lines', name: `BG: ${bgStem}`, line: { width: multiLineWidth, color: 'gray', dash: 'dash' }, opacity: 0.6 });
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

    const filesText = `${sorted.join(', ')} | BG: ${bgStem}`;
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

    const multiRawFilename = `${state.currentDate}_Multi_Raw_BG_${bgStem}`;
    Plotly.newPlot('plotMultiRaw', rawTraces, layoutRaw, getPlotConfig(multiRawFilename));
    document.getElementById('plotMultiRaw').dataset.filename = multiRawFilename;
    const multiDiffFilename = `${state.currentDate}_Multi_Diff_BG_${bgStem}`;
    Plotly.newPlot('plotMultiDiff', diffTraces, layoutDiff, getPlotConfig(multiDiffFilename));
    document.getElementById('plotMultiDiff').dataset.filename = multiDiffFilename;
  });

  /* ============================================================
     11. タブ5: 2範囲積算・比プロット
     ============================================================ */
  document.getElementById('btnPlotIntegration').addEventListener('click', async () => {
    const bgStem = document.getElementById('intBgSelect').value;
    if (!bgStem) { alert('Backgroundファイルを選択してください'); return; }
    const checked = Array.from(document.querySelectorAll('#intFileChecks input:checked')).map(c => c.value);
    if (checked.length === 0) { alert('積算対象の信号ファイルを1つ以上選択してください'); return; }

    const int1Start = parseFloat(document.getElementById('int1Start').value);
    const int1End = parseFloat(document.getElementById('int1End').value);
    const int2Start = parseFloat(document.getElementById('int2Start').value);
    const int2End = parseFloat(document.getElementById('int2End').value);
    const window_ = Math.max(1, parseInt(document.getElementById('intSmoothWin').value, 10) || 1);
    const xlabel = document.getElementById('intXLabel').value.trim() || 'X';

    const bgEntry = state.files[bgStem];
    const bgCsv = await ensureCsvParsed(bgStem);
    if (!bgCsv) { alert('BGのcsvファイルが見つかりません'); return; }
    applyGraphSize(['plotInt1', 'plotInt2', 'plotInt3', 'plotInt4']);
    const bgCond = bgEntry.condData || {};
    const bgTime = makeTimeAxis(bgCsv, bgCond);
    const bgCountsRaw = bgCsv.counts;
    const bgCountsSmooth = smooth(bgCountsRaw, window_);

    const results = [];
    for (const stem of checked) {
      const entry = state.files[stem];
      if (entry.labelType !== 'number') continue;
      const csv = await ensureCsvParsed(stem);
      if (!csv) continue;
      const cond = entry.condData || {};
      const sigTime = makeTimeAxis(csv, cond);
      const sigCountsRaw = csv.counts;
      const sigCountsSmooth = smooth(sigCountsRaw, window_);

      let diffCounts, bgRawInterp;
      if (arraysEqual(sigTime, bgTime)) {
        diffCounts = sigCountsSmooth.map((v, i) => v - bgCountsSmooth[i]);
        bgRawInterp = bgCountsRaw;
      } else {
        const bgSmoothInterp = interpLinear(sigTime, bgTime, bgCountsSmooth);
        diffCounts = sigCountsSmooth.map((v, i) => v - bgSmoothInterp[i]);
        bgRawInterp = interpLinear(sigTime, bgTime, bgCountsRaw);
      }

      const mask1 = maskRange(sigTime, int1Start, int1End);
      const mask2 = maskRange(sigTime, int2Start, int2End);

      const integral1 = sumMasked(diffCounts, mask1);
      const integral2 = sumMasked(diffCounts, mask2);
      const err1 = Math.sqrt(sumMasked(sigCountsRaw, mask1) + sumMasked(bgRawInterp, mask1));
      const err2 = Math.sqrt(sumMasked(sigCountsRaw, mask2) + sumMasked(bgRawInterp, mask2));

      results.push({ stem, label: entry.labelValue, integral1, integral2, err1, err2 });
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

    // 「⚙ 表示設定」のパレットが選択されていれば既定のblue/orange/greenをその配色に置き換え、
    // 線の太さも共通設定を反映する。
    const { palette: intPalette, lineWidth: intLineWidth } = getGraphStyle();
    function makeTrace(y, err, name, color, symbol, paletteIdx) {
      const resolvedColor = pickColor(intPalette, paletteIdx, color);
      const trace = { x: xs, y, type: 'scatter', mode, name, marker: { color: resolvedColor, symbol }, line: { color: resolvedColor, width: intLineWidth } };
      if (showErr) trace.error_y = { type: 'data', array: err, visible: true };
      return trace;
    }
    const filesText = `${results.map(r => r.stem).join(', ')} | BG: ${bgStem}`;
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

    const int1Filename = `${state.currentDate}_Integration_Range1_BG_${bgStem}`;
    Plotly.newPlot('plotInt1', [makeTrace(i1, e1, 'Range 1', 'blue', 'circle', 0)], baseLayout(`Range 1 (${int1Start} - ${int1End} ns)`, 'Integrated Counts'), getPlotConfig(int1Filename));
    document.getElementById('plotInt1').dataset.filename = int1Filename;
    const int2Filename = `${state.currentDate}_Integration_Range2_BG_${bgStem}`;
    Plotly.newPlot('plotInt2', [makeTrace(i2, e2, 'Range 2', 'orange', 'square', 0)], baseLayout(`Range 2 (${int2Start} - ${int2End} ns)`, 'Integrated Counts'), getPlotConfig(int2Filename));
    document.getElementById('plotInt2').dataset.filename = int2Filename;
    const int3Filename = `${state.currentDate}_Integration_Simultaneous_BG_${bgStem}`;
    Plotly.newPlot('plotInt3', [makeTrace(i1, e1, 'Range 1', 'blue', 'circle', 0), makeTrace(i2, e2, 'Range 2', 'orange', 'square', 1)], baseLayout('Simultaneous Plot', 'Integrated Counts'), getPlotConfig(int3Filename));
    document.getElementById('plotInt3').dataset.filename = int3Filename;
    const int4Filename = `${state.currentDate}_Integration_Ratio_BG_${bgStem}`;
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

  /* ============================================================
     12. 「5. コードを編集する」— app.js の表示・編集・プレビュー再実行
     ============================================================
     index.html と app.js を fetch で取得し、テキストエリアに表示する。
     「変更を反映して再実行」を押すと、取得済みの index.html の <script src="app.js">
     をテキストエリアの内容に差し替えたHTMLを組み立て、iframeのsrcdocに設定して
     独立したJS実行環境（別のグローバルスコープ）でツール全体を再読み込みする。
     こうすることで、このページ自体の const/関数宣言と衝突せずに編集結果をその場で試せる。
     fetch は同一オリジンの静的ファイル読込みが前提のため、file:// で直接開いた場合は
     失敗することがある（その場合はローカルサーバー経由か公開ページの利用を案内する）。 */
  (function setupCodeEditor() {
    const codeEditor = document.getElementById('codeEditor');
    const codeStatus = document.getElementById('codeStatus');
    const previewFrame = document.getElementById('codePreviewFrame');
    const btnRunCode = document.getElementById('btnRunCode');
    const btnResetCode = document.getElementById('btnResetCode');
    const btnDownloadCode = document.getElementById('btnDownloadCode');
    if (!codeEditor || !previewFrame || !btnRunCode) return;

    let originalCode = null;
    let originalIndexHtml = null;

    function setStatus(msg, isWarn) {
      codeStatus.textContent = msg || '';
      codeStatus.classList.toggle('warn', !!isWarn);
    }

    async function loadOriginalSource() {
      setStatus('読み込み中...');
      try {
        const [codeText, htmlText] = await Promise.all([
          fetch('app.js', { cache: 'no-store' }).then(r => { if (!r.ok) throw new Error(String(r.status)); return r.text(); }),
          fetch('index.html', { cache: 'no-store' }).then(r => { if (!r.ok) throw new Error(String(r.status)); return r.text(); }),
        ]);
        originalCode = codeText;
        originalIndexHtml = htmlText;
        codeEditor.value = originalCode;
        setStatus('');
      } catch (e) {
        setStatus('コードの読み込みに失敗しました。file:// で直接開いている場合はローカルサーバー経由で開くか、公開ページをご利用ください。', true);
      }
    }

    // 取得済みの index.html から「5. コードを編集する」欄自体を取り除き、
    // <script src="app.js"> を編集後のコードのインラインscriptに差し替えたHTML文字列を作る。
    function buildPreviewHtml(code) {
      const doc = new DOMParser().parseFromString(originalIndexHtml, 'text/html');
      const codeSection = doc.getElementById('section-code');
      if (codeSection) codeSection.remove();
      doc.querySelectorAll('script[src="app.js"]').forEach(scriptEl => {
        const inline = doc.createElement('script');
        inline.textContent = code;
        scriptEl.replaceWith(inline);
      });
      return '<!DOCTYPE html>\n' + doc.documentElement.outerHTML;
    }

    btnRunCode.addEventListener('click', () => {
      if (originalIndexHtml == null) { alert('元のコードがまだ読み込めていません。'); return; }
      try {
        previewFrame.srcdoc = buildPreviewHtml(codeEditor.value);
        setStatus(`反映しました（${new Date().toLocaleTimeString()}）`);
      } catch (e) {
        setStatus('プレビューの生成に失敗しました: ' + (e && e.message ? e.message : e), true);
      }
    });

    btnResetCode.addEventListener('click', () => {
      if (originalCode == null) { alert('元のコードがまだ読み込めていません。'); return; }
      if (codeEditor.value !== originalCode && !confirm('編集内容を破棄して元のコードに戻しますか？')) return;
      codeEditor.value = originalCode;
      setStatus('元のコードに戻しました');
    });

    btnDownloadCode.addEventListener('click', () => {
      const blob = new Blob([codeEditor.value], { type: 'text/javascript' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'app.js';
      document.body.appendChild(a);
      a.click();
      a.remove();
      URL.revokeObjectURL(url);
    });

    loadOriginalSource();
  })();

  /* ---------- 初期描画（IndexedDBに保存済みのファイル・設定があれば自動復元する） ---------- */
  renderAll();
  restoreFromDB();
}
