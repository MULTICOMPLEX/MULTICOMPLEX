"use strict";/*Compiled using Cheerp (R) by Leaning Technologies Ltd*/var x=Math.imul;var y=Math.fround;var oSlot=0;var nullArray=[null];var nullObj={d:nullArray,o:0};function v(){var b=null,d=null,a=null;a="cheerpWorker.js";b=new Worker(a);a="message";d=u;b.addEventListener(a,d);a="Hello World";b.postMessage(a);}function u(a){var b=null;b=a.data;console.log(b);}function e(f,g){var a=null,d=0,b=null;a=String();if((f[g]&255)===0)return String(a);d=0;while(1){b=String.fromCharCode(f[g+d|0]<<24>>24);a=a.concat(b);d=d+1|0;if((f[g+d|0]&255)!==0)continue;break;}return String(a);}var j=new Uint8Array([99,104,101,101,114,112,87,111,114,107,101,114,46,106,115,0]);var k=new Uint8Array([109,101,115,115,97,103,101,0]);var i=new Uint8Array([72,101,108,108,111,32,87,111,114,108,100,0]);v();