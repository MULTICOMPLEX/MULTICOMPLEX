"use strict";/*Compiled using Cheerp (R) by Leaning Technologies Ltd*/var uX=Math.imul;var uY=Math.fround;var oSlot=0;var nullArray=[null];var nullObj={d:nullArray,o:0};function vm(p){var b=null,f='function';if(typeof fetch===f)b=fetch(p).then(r=>r.arrayBuffer());else if(typeof require===f){p=require('path').join(__dirname, p);b=new Promise((y,n)=>{require('fs').readFile(p,(e,d)=>{if(e)n(e);else y(d);});});}else b=new Promise((y,n)=>{y(read(p,'binary'));});return b;}function g$(){return +Date.now();}function g8(){return ~~( +Math.random()*4294967296)|0;}function ml(l,m,j){var i=0,g=null;g=lM(l,m,j);i=j-1|0;if((l[m+i|0]&255)===10){g=g.substr(0,i);console.log(g);return;}console.log(g);}function lM(t,u,q){var l=0,n=0,g=null,o=0,p=0,i=0,j=null;g=String();if((q|0)===0)return g;p=q;o=0;while(1){l=t[u+o|0]|0;if((l&255)!==0){n=l&255;if(l<<24>-16777216){i=n;}else if((l&255)<192){i=n&63|i<<6;}else if((l&255)<224){i=n&31;}else if((l&255)<240){i=n&15;}else{i=n&7;}p=p-1|0;o=o+1|0;a:{if((p|0)!==0)if((t[u+o|0]&192)===128)break a;if(i>>>0<65536){j=String.fromCharCode(i);g=g.concat(j);}else{i=i-65536|0;j=String.fromCharCode((i>>>10)+55296|0);g=g.concat(j);j=String.fromCharCode((i&1023)+56320|0);g=g.concat(j);}}if((p|0)!==0)continue;return g;}break;}return g;}function mm(){throw new Error("Abort called");}function uu(l,j){var i=0,g=0;if((j|0)<2)return 1|0;g=1;i=j;while(1){g=uX(g,i)|0;if((i|0)<3)return g|0;i=i-1|0;continue;}}function uv(g){return g.a2;}function ut(n,l){var j=0,i=-0.,g=null;n.d0=l;n.i1=0;j=1052936>>0;g=a;g=g;i=+g.BYTES_PER_ELEMENT;n.a2=new Float64Array(g.buffer,(+(uX(j,~~i)>>>0)),4);}function uy(){var g=null,i=null;g=dO(a,(1069991>>0));i=ll;document.addEventListener(g,i);}function dO(g,h){var i=null,l=0,j=null;i=String();if((g[h]&255)===0)return String(i);l=0;while(1){j=String.fromCharCode(g[h+l|0]<<24>>24);i=i.concat(j);l=l+1|0;if((g[h+l|0]&255)!==0)continue;break;}return String(i);}function lk(){var n=null,l=null,q=null,p=-0.,o=-0.,i=null,j=null,k=0,g=0;n=fq();l=-64+n|0;d3(l);document.body;if(!(iE|0)){i=dO(a,(1070634>>0));j=document.getElementById(i);iD=j;iE=1;}i=dO(a,(1070573>>0));q=document.getElementById(i);i=dO(a,(1070742>>0));j=document.getElementById(i);i=dO(a,(1070751>>0));document.getElementById(i);i=dO(a,(1070757>>0));document.getElementById(i);i=q.value;j=j.value;p=+parseInt(i);o=+parseInt(j);f[1048856>>3]=p/1000;f[1048864>>3]=o/1000;hz(hy(hz(1077880|0,1070763|0)|0|0,1048840|0)|0|0,1070794|0)|0;hx();g=32+l|0;fX(g,1048840|0);hy(1077880|0,g)|0;i=iD;g=16+l|0;lh(g);k=c[8+g>>2];j=a;j=dO(j,k);i.textContent=j;g7(g);g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);lf(g);g=l|0;le(g);ld(g);g7(g);d3(n);}function lm(){var i=null,g=0;g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);c[8+g>>2]=12;i=gc;+requestAnimationFrame(i);}function gc(){var g=null;lk();g=gc;+requestAnimationFrame(g);}function ll(){var i=null,g=0;g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1077872>>2]|0)>>2]|0;g=(1077872+g|0);c[8+g>>2]=12;i=gc;+requestAnimationFrame(i);}function uz(D,I){var y=null,j=0,x=0,o=0,p=0,q=0,B=0,z=-0.,w=0,i=null,l=null,g=0,t=-0.,v=0,n=0;y=fq();i=-57736+y|0;d3(i);j=57704+i|0;fR(j);x=57696+i|0;la(x);o=96+i|0;if(!(ic|0)){g=o>>0;l=a;l=l;t=+l.BYTES_PER_ELEMENT;ib=new Float64Array(l.buffer,(+(uX(g,~~t)>>>0)),7200);ic=1;}l=o|0;g=64+i|0;p=32+i|0;q=i|0;v=0;t=-1;while(1){B=uX(v,200)|0;n=0;while(1){z=(+(n-100|0))/100;f[16+j>>3]=z;f[24+j>>3]=t;fR(g);hp(p,j,I);ho(g,p);w=n+B|0;f[(w<<3)+l>>3]=+f[16+g>>3];f[(w<<3)+28800+l>>3]=+f[24+g>>3];f[16+j>>3]=t;f[24+j>>3]=z;hp(q,j,I);ho(g,q);f[(w<<3)+14400+l>>3]=+f[16+g>>3];f[(w<<3)+43200+l>>3]=+f[24+g>>3];n=n+1|0;if((n|0)!==200)continue;break;}v=v+1|0;if((v|0)===9){i=ib;d3(y);return i;}t+=.25;continue;}}function ux(n,l,j,i,g){return (l-j)/(i-j);}function uA(t,q,p){var n=null,i=0,o=-0.,j=null,g=0,l=null;n=fq();j=-64+n|0;d3(j);if(!(ia|0)){g=1052936>>0;l=a;l=l;o=+l.BYTES_PER_ELEMENT;h$=new Float64Array(l.buffer,(+(uX(g,~~o)>>>0)),4);ia=1;}g=32+j|0;fR(g);f[16+g>>3]=q;f[24+g>>3]=p;i=j|0;fX(i,g);f[1052936>>3]=+f[16+i>>3];f[1052944>>3]=+f[24+i>>3];j=h$;d3(n);return j;}function uw(g){}function i9(i,g){i=i|0;g=g|0;ml(a,(i>>0),g);}var iE=0;var iD=null;var ic=0;var ib=null;var ia=0;var h$=null;function um(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a4;obj.a4.o=2;obj.a4.a=a;return obj;}function ul(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a3;obj.a3.o=2;obj.a3.a=a;return obj;}function vh(bytes){var pages=(bytes+65535)>>16;try{__asm.vk.grow(pages);vi(__asm.vk.buffer);return pages<<16;}catch(e){return -1;}}var a=null,c=null,f=null,__asm=null,__heap=null;function vj(){throw new Error('this should be unreachable');};var hz=null;var hy=null;var hx=null;var fX=null;var lh=null;var g7=null;var lf=null;var le=null;var ld=null;var fR=null;var la=null;var hp=null;var ho=null;var fq=null;var d3=null;var Graphics=vj;var JsStruct=vj;vj.promise=vm('derivatives.wasm').then(g=>WebAssembly.instantiate(g,{i:{a3:vj,a2:vj,kK:vj,g8:g8,g$:g$,lm:lm,mm:mm,i9:i9,u9:Math.pow,vh:vh,}})).then(g=>{__asm=g.instance.exports;__heap=__asm.vk.buffer;vi(__heap);hz=__asm.hz;hy=__asm.hy;hx=__asm.hx;fX=__asm.fX;lh=__asm.lh;g7=__asm.g7;lf=__asm.lf;le=__asm.le;ld=__asm.ld;fR=__asm.fR;la=__asm.la;hp=__asm.hp;ho=__asm.ho;fq=__asm.fq;d3=__asm.d3;Graphics=function (){this.i0=0;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}uw(this);};Graphics.prototype.update_array=function (a0,a1){return uA(this,a0,a1);};Graphics.prototype.normalize=function (a0,a1,a2,a3){return ux(this,a0,a1,a2,a3);};Graphics.prototype.conformal_map=function (a0){return uz(this,a0);};Graphics.loadCallback=function (){return ll();};Graphics.initialize=function (){return lm();};Graphics.mainLoop=function (){return lk();};Graphics.init=function (){return uy();};JsStruct=function (a0){this.d0=uY(0.);this.i1=0;this.a2=null;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}ut(this,a0);};JsStruct.prototype.test=function (){return uv(this);};JsStruct.prototype.factorial=function (a0){return uu(this,a0);};Graphics.promise=JsStruct.promise=Promise.resolve();__asm.h9();__asm.h8();__asm.kY();});function vi(g){a=new Uint8Array(g);c=new Int32Array(g);f=new Float64Array(g);}