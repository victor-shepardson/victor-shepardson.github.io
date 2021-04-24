// TODO: run shader on canvas
// TODO: feedback shader
// TODO: fix blank bar at top
// TODO: use canvas as background / float text over canvas?
// TODO: avoid page load?

var rasterizeHTML = require('rasterizehtml');
var fs = require('fs');

var canvas = document.getElementById("canvas"),
  // context = canvas.getContext('2d'),
  html = fs.readFileSync(__dirname + '/index.html', 'utf8');
//
// rasterizeHTML.drawHTML(html).then(function(renderResult) {
//   // canvas.width = window.innerWidth;
//   // canvas.height = window.innerHeight;
//   context.drawImage(renderResult.image, 0, 0);
// });

canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
rasterizeHTML.drawHTML(html, canvas);