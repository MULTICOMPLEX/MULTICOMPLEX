"use strict";/*Compiled using Cheerp (R) by Leaning Technologies Ltd*/var s5=Math.imul;var s6=Math.fround;var oSlot=0;var nullArray=[null];var nullObj={d:nullArray,o:0};function tu(p){var b=null,f='function';if(typeof fetch===f)b=fetch(p).then(r=>r.arrayBuffer());else if(typeof require===f){p=require('path').join(__dirname, p);b=new Promise((y,n)=>{require('fs').readFile(p,(e,d)=>{if(e)n(e);else y(d);});});}else b=new Promise((y,n)=>{y(read(p,'binary'));});return b;}function fp(){return +Date.now();}function fm(){return ~~( +Math.random()*4294967296)|0;}function kp(l,m,j){var i=0,g=null;g=ko(l,m,j);i=j-1|0;if((l[m+i|0]&255)===10){g=g.substr(0,i);console.log(g);return;}console.log(g);}function ko(r,s,q){var l=0,n=0,g=null,o=0,p=0,i=0,j=null;g=String();if((q|0)===0)return g;p=q;o=0;while(1){l=r[s+o|0]|0;if((l&255)!==0){n=l&255;if(l<<24>-16777216){i=n;}else if((l&255)<192){i=n&63|i<<6;}else if((l&255)<224){i=n&31;}else if((l&255)<240){i=n&15;}else{i=n&7;}p=p-1|0;o=o+1|0;a:{if((p|0)!==0)if((r[s+o|0]&192)===128)break a;if(i>>>0<65536){j=String.fromCharCode(i);g=g.concat(j);}else{i=i-65536|0;j=String.fromCharCode((i>>>10)+55296|0);g=g.concat(j);j=String.fromCharCode((i&1023)+56320|0);g=g.concat(j);}}if((p|0)!==0)continue;return g;}break;}return g;}function kq(){throw new Error("Abort called");}function sE(l,j){var i=0,g=0;if((j|0)<2)return 1|0;g=1;i=j;while(1){g=s5(g,i)|0;if((i|0)<3)return g|0;i=i-1|0;continue;}}function sF(g){return g.a2;}function sD(n,l){var j=0,i=-0.,g=null;n.d0=l;n.i1=0;j=1052936>>0;g=a;g=g;i=+g.BYTES_PER_ELEMENT;n.a2=new Float64Array(g.buffer,(+(s5(j,~~i)>>>0)),4);}function sH(){var g=null,i=null;g=c$(a,(1070007>>0));i=jV;document.addEventListener(g,i);}function c$(g,h){var i=null,l=0,j=null;i=String();if((g[h]&255)===0)return String(i);l=0;while(1){j=String.fromCharCode(g[h+l|0]<<24>>24);i=i.concat(j);l=l+1|0;if((g[h+l|0]&255)!==0)continue;break;}return String(i);}function jU(){var p=null,l=null,t=null,r=-0.,q=-0.,o=0,n=0,i=null,j=null,k=0,g=0;p=fe();l=-3104+p|0;d1(l);he=document.body;if(!(hd|0)){i=c$(a,(1070650>>0));j=document.createElement(i);e8=j;hd=1;}i=c$(a,(1070589>>0));t=document.getElementById(i);i=c$(a,(1070748>>0));j=document.getElementById(i);i=c$(a,(1070757>>0));document.getElementById(i);i=c$(a,(1070763>>0));document.getElementById(i);i=t.value;j=j.value;r=+parseInt(i);q=+parseInt(j);f[1048856>>3]=r/100;f[1048864>>3]=q/100;gf(ge(gf(1077880|0,1070769|0)|0|0,1048840|0)|0|0,1070793|0)|0;eV();g=3072+l|0;o=1552+l|0;n=32+l|0;eU(n);gd(o,n);eT(g,o);jR(g);ge(1077880|0,1053416|0)|0;i=e8;g=16+l|0;jQ(g);k=c[8+g>>2];j=a;j=c$(j,k);i.textContent=j;fl(g);g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);jO(g);g=l|0;jN(g);jM(g);fl(g);he.appendChild(e8);d1(p);}function jW(){var i=null,g=0;g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);c[8+g>>2]=12;i=eQ;+requestAnimationFrame(i);}function eQ(){var g=null;jU();g=eQ;+requestAnimationFrame(g);}function jV(){var i=null,g=0;g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);c[8+g>>2]=12;i=eQ;+requestAnimationFrame(i);}function sI(q){var o=null,l=0,j=0,p=-0.,i=null,g=0,n=null;o=fe();i=-3072+o|0;d1(i);if(!(gT|0)){g=1052936>>0;n=a;n=n;p=+n.BYTES_PER_ELEMENT;gS=new Float64Array(n.buffer,(+(s5(g,~~p)>>>0)),4);gT=1;}eV();g=3040+i|0;l=1520+i|0;j=i|0;eU(j);gd(l,j);eT(g,l);f[1052936>>3]=+f[16+g>>3];f[1052944>>3]=+f[24+g>>3];i=gS;d1(o);return i;}function sG(g){}function hJ(i,g){i=i|0;g=g|0;kp(a,(i>>0),g);}var he=null;var hd=0;var e8=null;var gT=0;var gS=null;function sw(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a4;obj.a4.o=2;obj.a4.a=a;return obj;}function sv(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a3;obj.a3.o=2;obj.a3.a=a;return obj;}function tp(bytes){var pages=(bytes+65535)>>16;try{__asm.ts.grow(pages);tq(__asm.ts.buffer);return pages<<16;}catch(e){return -1;}}var a=null,c=null,f=null,__asm=null,__heap=null;function tr(){throw new Error('this should be unreachable');};var gf=null;var ge=null;var eV=null;var eU=null;var gd=null;var eT=null;var jR=null;var jQ=null;var fl=null;var jO=null;var jN=null;var jM=null;var fe=null;var d1=null;var Graphics=tr;var JsStruct=tr;tr.promise=tu('derivatives.wasm').then(g=>WebAssembly.instantiate(g,{i:{aB:tr,aA:tr,kk:tr,jW:jW,fm:fm,fp:fp,kq:kq,hJ:hJ,tp:tp,}})).then(g=>{__asm=g.instance.exports;__heap=__asm.ts.buffer;tq(__heap);gf=__asm.gf;ge=__asm.ge;eV=__asm.eV;eU=__asm.eU;gd=__asm.gd;eT=__asm.eT;jR=__asm.jR;jQ=__asm.jQ;fl=__asm.fl;jO=__asm.jO;jN=__asm.jN;jM=__asm.jM;fe=__asm.fe;d1=__asm.d1;Graphics=function (){this.i0=0;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}sG(this);};Graphics.prototype.update_array=function (){return sI(this);};Graphics.loadCallback=function (){return jV();};Graphics.initialize=function (){return jW();};Graphics.mainLoop=function (){return jU();};Graphics.init=function (){return sH();};JsStruct=function (a0){this.d0=s6(0.);this.i1=0;this.a2=null;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}sD(this,a0);};JsStruct.prototype.test=function (){return sF(this);};JsStruct.prototype.factorial=function (a0){return sE(this,a0);};Graphics.promise=JsStruct.promise=Promise.resolve();__asm.gQ();__asm.gP();__asm.jy();});function tq(g){a=new Uint8Array(g);c=new Int32Array(g);f=new Float64Array(g);}