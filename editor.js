/* ============================================================
   コード編集ページ（editor.html）専用のスクリプト。
   index.html と app.js を fetch で取得し、テキストエリアに app.js の内容を表示する。
   「変更を反映して再実行」を押すと、取得済みの index.html の <script src="app.js">
   をテキストエリアの内容に差し替えたHTMLを組み立て、iframeのsrcdocに設定して
   独立したJS実行環境（別のグローバルスコープ）でツール全体を再読み込みする。
   fetch は同一オリジンの静的ファイル読込みが前提のため、file:// で直接開いた場合は
   失敗することがある（その場合はローカルサーバー経由か公開ページの利用を案内する）。
   ============================================================ */
(function setupCodeEditor() {
  const codeEditor = document.getElementById('codeEditor');
  const codeStatus = document.getElementById('codeStatus');
  const previewFrame = document.getElementById('codePreviewFrame');
  const btnRunCode = document.getElementById('btnRunCode');
  const btnResetCode = document.getElementById('btnResetCode');
  const btnDownloadCode = document.getElementById('btnDownloadCode');

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

  // 取得済みの index.html の <script src="app.js"> を、編集後のコードのインラインscriptに差し替えたHTML文字列を作る。
  function buildPreviewHtml(code) {
    const doc = new DOMParser().parseFromString(originalIndexHtml, 'text/html');
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
