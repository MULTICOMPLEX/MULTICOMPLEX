"use strict";/*Compiled using Cheerp (R) by Leaning Technologies Ltd*/var pU=Math.imul;var pV=Math.fround;var oSlot=0;var nullArray=[null];var nullObj={d:nullArray,o:0};function qj(p){var b=null,f='function';if(typeof fetch===f)b=fetch(p).then(r=>r.arrayBuffer());else if(typeof require===f){p=require('path').join(__dirname, p);b=new Promise((y,n)=>{require('fs').readFile(p,(e,d)=>{if(e)n(e);else y(d);});});}else b=new Promise((y,n)=>{y(read(p,'binary'));});return b;}function hi(){return +Date.now();}function hF(q,r,p){var f=null,m=0,o=0,n=0,j=0,h=0,l=null;f=String();a:if((p|0)!==0){o=p;m=0;while(1){j=q[r+m|0]|0;if((j&255)!==0){while(1){h=j&255;if(j<<24<=-16777216)if((j&255)<192){h=h&63|n<<6;}else if((j&255)<224){h&=31;}else if((j&255)<240){h&=15;}else{h&=7;}o=o-1|0;m=m+1|0;if((o|0)!==0){j=q[r+m|0]|0;if((j&192)===128){if((j&255)===0)break a;n=h;continue;}j=0;}else{j=1;o=0;}break;}if(h>>>0<65536){n=h;}else{n=h-65536|0;l=String.fromCharCode((n>>>10)+55296|0);f=f.concat(l);h=(n&1023)+56320|0;}l=String.fromCharCode(h);f=f.concat(l);if(!(j))continue;}break;}}m=p-1|0;if((q[r+m|0]&255)===10){l=f.substr(0,m);console.log(l);return;}console.log(f);}function hO(){throw new Error("Abort called");}function i2(){var h=null,f=0;es();f=c[16+(c[1074104>>2]|0)>>2]|0;f=(1074104+f|0);c[4+f>>2]=c[4+f>>2]& -261|4;f=c[16+(c[1074104>>2]|0)>>2]|0;f=(1074104+f|0);c[8+f>>2]=12;h=er;+requestAnimationFrame(h);}function er(){var f=null;i0();f=er;+requestAnimationFrame(f);}function i0(){var o=null,l=null,s=null,q=-0.,p=-0.,n=0,m=0,h=null,i=0,j=null,k=0,f=0;o=gd();l=-3104+o|0;dO(l);if(!(eJ|0)){eK=document.body;eJ=1;}if(!(eL|0)){h="h1";j=document.createElement(h);dG=j;eL=1;}h="myRange1";s=document.getElementById(h);h="myRange2";j=document.getElementById(h);h="demo1";document.getElementById(h);h="demo2";document.getElementById(h);h=s.value;j=j.value;q=+parseInt(h);p=+parseInt(j);e[16+1054104>>3]=q/100;e[24+1054104>>3]=p/100;et(en(et(8+1074104|0,1071452|0)|0|0,1054104|0)|0|0,1070877|0)|0;es();n=3072+l|0;f=1552+l|0;m=32+l|0;eI(m);iZ(f,m);eE(n,f);en(8+1074104|0,n)|0;h=dG;f=16+l|0;iX(f);k=c[8+f>>2];j=a;j=cu(j,k);h.textContent=j;d0(f);f=c[16+(c[1074104>>2]|0)>>2]|0;f=(1074104+f|0);i=c[24+f>>2];h=a;c[16+f>>2]=h===nullArray&&i===0?1:0;f=l|0;hj(f);iW(f);d0(f);eK.appendChild(dG);dO(o);}function cu(f,g){var j=null,h=0,m=0,l=null;j=String();h=f[g]|0;if((h&255)===0)return String(j);m=0;while(1){l=String.fromCharCode(h<<24>>24);j=j.concat(l);m=m+1|0;h=f[g+m|0]|0;if((h&255)!==0)continue;break;}return String(j);}function gc(h,f){h=h|0;f=f|0;hF(a,(h>>0),f);}var eJ=0;var eK=null;var eL=0;var dG=null;function pw(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a4;obj.a4.o=2;obj.a4.a=a;return obj;}function pv(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a3;obj.a3.o=2;obj.a3.a=a;return obj;}function qe(bytes){var pages=(bytes+65535)>>16;try{__asm.qh.grow(pages);qf(__asm.qh.buffer);return pages<<16;}catch(e){return -1;}}function qf(f){a=new Uint8Array(f);b=new Uint16Array(f);c=new Int32Array(f);d=new Float32Array(f);e=new Float64Array(f);}var a=null,b=null,c=null,d=null,e=null,__asm=null,__heap=null;function qg(){throw new Error('this should be unreachable');};var es=null;var et=null;var en=null;var eI=null;var iZ=null;var eE=null;var iX=null;var d0=null;var hj=null;var iW=null;var gd=null;var dO=null;qj('derivatives.wasm').then(f=>WebAssembly.instantiate(f,{i:{i2:i2,hi:hi,hO:hO,gc:gc,qe:qe,}}),console.log).then(f=>{__asm=f.instance.exports;__heap=__asm.qh.buffer;qf(__heap);es=__asm.es;et=__asm.et;en=__asm.en;eI=__asm.eI;iZ=__asm.iZ;eE=__asm.eE;iX=__asm.iX;d0=__asm.d0;hj=__asm.hj;iW=__asm.iW;gd=__asm.gd;dO=__asm.dO;Promise.resolve();__asm.fa();__asm.i3();},console.log,console.log);