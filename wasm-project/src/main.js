fetch('../out/main.wasm').then(response =>
  response.arrayBuffer()
).then(bytes => WebAssembly.instantiate(bytes)).then(results => {
  instance = results.instance;
  
  
  document.getElementById("add").innerHTML = "1 + 2 = " + instance.exports.add(1,2) + "<br />";
  document.getElementById("sub").innerHTML = "1 - 2 = " + instance.exports.sub(1,2);

}).catch(console.error);
