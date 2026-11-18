// If the URL fragment points inside a collapsed <details> element (e.g., a
// link to a deprecated metadata field), open it and scroll to the target.
function openDetailsForHash() {
  const id = decodeURIComponent(location.hash.slice(1));
  if (!id) return;
  const target = document.getElementById(id);
  if (!target) return;
  let details = target.closest("details");
  if (!details || details.open) return;
  while (details) {
    details.open = true;
    details = details.parentElement.closest("details");
  }
  target.scrollIntoView();
}

window.addEventListener("hashchange", openDetailsForHash);
document$.subscribe(openDetailsForHash);
