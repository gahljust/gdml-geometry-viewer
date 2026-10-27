/* One explicit JSON download; no background saving or persistent cache. */
(() => {
  const button = document.getElementById('save-view');
  button.onclick = () => {
    const detail = {snapshot: null, error: ''};
    document.dispatchEvent(new CustomEvent('save-current-view', {detail}));
    if (!detail.snapshot) {
      button.textContent = detail.error || 'Load a view first';
      setTimeout(() => { button.textContent = 'Save view'; }, 2000);
      return;
    }
    const data = {schema: 'moller-view-v1', saved_utc: new Date().toISOString(), ...detail.snapshot};
    const url = URL.createObjectURL(new Blob([JSON.stringify(data)], {type: 'application/json'}));
    const link = document.createElement('a');
    link.href = url;
    link.download = `moller-${data.page}-view.json`;
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  };
})();
