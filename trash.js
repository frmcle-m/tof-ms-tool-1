/* ============================================================
   ゴミ箱ページ（trash.html）専用のスクリプト。
   解析ツール本体（app.js）と同じ IndexedDB（tof_tool_db）を直接読み書きする。
   DB_NAME・DB_VERSION・ストア名は app.js 側と必ず一致させること。
   ============================================================ */

// calfact は整数値のため、小数点以下は切り捨てて表示する（app.js と同じロジック）
function formatCalfact(v) {
  if (v === undefined || v === null || v === '') return v;
  const n = parseFloat(v);
  if (!Number.isFinite(n)) return v;
  return String(Math.trunc(n));
}

const DB_NAME = 'tof_tool_db';
const DB_VERSION = 4;
const STORE_FILES = 'files';
const STORE_META = 'meta';
const STORE_TRASH = 'trash';

function openDB() {
  return new Promise((resolve) => {
    if (typeof indexedDB === 'undefined') { resolve(null); return; }
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

function formatTrashedAt(ts) {
  if (!ts) return '';
  try {
    return new Date(ts).toLocaleString('ja-JP');
  } catch (e) {
    return '';
  }
}

function describeLabelType(labelType) {
  if (labelType === 'number') return '数値';
  if (labelType === 'string') return '文字列';
  return '(未設定)';
}

async function restoreRecord(stem) {
  const records = await dbGetAll(STORE_TRASH);
  const record = records.find(r => r.stem === stem);
  if (!record) return;
  const { trashedAt, ...fileRecord } = record;
  await dbPut(STORE_FILES, fileRecord);
  await dbDeleteKey(STORE_TRASH, stem);
  await render();
}

async function permanentlyDelete(stem) {
  if (!confirm(`「${stem}」を完全に削除します。この操作は元に戻せません。続けますか？`)) return;
  await dbDeleteKey(STORE_TRASH, stem);
  await render();
}

async function emptyTrash() {
  const records = await dbGetAll(STORE_TRASH);
  if (records.length === 0) { alert('ゴミ箱は空です。'); return; }
  if (!confirm(`ゴミ箱の中の${records.length}件のファイルをすべて完全に削除します。この操作は元に戻せません。続けますか？`)) return;
  for (const r of records) {
    await dbDeleteKey(STORE_TRASH, r.stem);
  }
  await render();
}

async function render() {
  const records = await dbGetAll(STORE_TRASH);
  records.sort((a, b) => (b.trashedAt || 0) - (a.trashedAt || 0));

  const table = document.getElementById('trashTable');
  const empty = document.getElementById('trashEmpty');
  const countText = document.getElementById('trashCountText');
  const tbody = document.getElementById('trashTableBody');

  countText.textContent = records.length > 0 ? `${records.length}件のファイルがゴミ箱に入っています。` : '';

  if (records.length === 0) {
    table.classList.add('hidden');
    empty.classList.remove('hidden');
    return;
  }
  table.classList.remove('hidden');
  empty.classList.add('hidden');
  tbody.innerHTML = '';

  for (const r of records) {
    const cond = r.condData || {};
    const tr = document.createElement('tr');

    const tdName = document.createElement('td');
    tdName.textContent = r.stem;
    const tdDate = document.createElement('td');
    tdDate.textContent = r.dateKey || '';
    const tdLabel = document.createElement('td');
    tdLabel.textContent = r.label || '';
    const tdType = document.createElement('td');
    tdType.textContent = describeLabelType(r.labelType);
    const tdCond = document.createElement('td');
    tdCond.textContent = r.hasCond ? '✓' : '✗';
    const tdCsv = document.createElement('td');
    tdCsv.textContent = r.hasCsv ? '✓' : '✗';
    const tdTime = document.createElement('td');
    tdTime.textContent = cond.measurement_time || '';
    const tdSweeps = document.createElement('td');
    tdSweeps.textContent = cond.SWEEPS || '';
    const tdCalfact = document.createElement('td');
    tdCalfact.textContent = ('calfact' in cond) ? formatCalfact(cond.calfact) : '-';
    const tdTrashedAt = document.createElement('td');
    tdTrashedAt.textContent = formatTrashedAt(r.trashedAt);

    const tdActions = document.createElement('td');
    const btnRestore = document.createElement('button');
    btnRestore.type = 'button';
    btnRestore.className = 'btn btn-ghost btn-small';
    btnRestore.textContent = '元に戻す';
    btnRestore.addEventListener('click', () => restoreRecord(r.stem));
    const btnDelete = document.createElement('button');
    btnDelete.type = 'button';
    btnDelete.className = 'btn btn-ghost btn-small btn-danger';
    btnDelete.textContent = '完全に削除';
    btnDelete.addEventListener('click', () => permanentlyDelete(r.stem));
    tdActions.appendChild(btnRestore);
    tdActions.appendChild(btnDelete);

    tr.appendChild(tdName); tr.appendChild(tdDate); tr.appendChild(tdLabel); tr.appendChild(tdType);
    tr.appendChild(tdCond); tr.appendChild(tdCsv); tr.appendChild(tdTime); tr.appendChild(tdSweeps);
    tr.appendChild(tdCalfact); tr.appendChild(tdTrashedAt); tr.appendChild(tdActions);
    tbody.appendChild(tr);
  }
}

document.getElementById('btnEmptyTrash').addEventListener('click', emptyTrash);

render();
