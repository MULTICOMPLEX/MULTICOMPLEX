"use strict";/*Compiled using Cheerp (R) by Leaning Technologies Ltd*/var uZ=Math.imul;var u0=Math.fround;var oSlot=0;var nullArray=[null];var nullObj={d:nullArray,o:0};function vo(p){var b=null,f='function';if(typeof fetch===f)b=fetch(p).then(r=>r.arrayBuffer());else if(typeof require===f){p=require('path').join(__dirname, p);b=new Promise((y,n)=>{require('fs').readFile(p,(e,d)=>{if(e)n(e);else y(d);});});}else b=new Promise((y,n)=>{y(read(p,'binary'));});return b;}function g8(){return +Date.now();}function g5(){return ~~( +Math.random()*4294967296)|0;}function k4(l,m,j){var i=0,g=null;g=k3(l,m,j);i=j-1|0;if((l[m+i|0]&255)===10){g=g.substr(0,i);console.log(g);return;}console.log(g);}function k3(t,u,q){var l=0,n=0,g=null,o=0,p=0,i=0,j=null;g=String();if((q|0)===0)return g;p=q;o=0;while(1){l=t[u+o|0]|0;if((l&255)!==0){n=l&255;if(l<<24>-16777216){i=n;}else if((l&255)<192){i=n&63|i<<6;}else if((l&255)<224){i=n&31;}else if((l&255)<240){i=n&15;}else{i=n&7;}p=p-1|0;o=o+1|0;a:{if((p|0)!==0)if((t[u+o|0]&192)===128)break a;if(i>>>0<65536){j=String.fromCharCode(i);g=g.concat(j);}else{i=i-65536|0;j=String.fromCharCode((i>>>10)+55296|0);g=g.concat(j);j=String.fromCharCode((i&1023)+56320|0);g=g.concat(j);}}if((p|0)!==0)continue;return g;}break;}return g;}function k5(){throw new Error("Abort called");}function uw(l,j){var i=0,g=0;if((j|0)<2)return 1|0;g=1;i=j;while(1){g=uZ(g,i)|0;if((i|0)<3)return g|0;i=i-1|0;continue;}}function ux(g){return g.a2;}function uv(n,l){var j=0,i=-0.,g=null;n.d0=l;n.i1=0;j=1053104>>0;g=a;g=g;i=+g.BYTES_PER_ELEMENT;n.a2=new Float64Array(g.buffer,(+(uZ(j,~~i)>>>0)),4);}function uA(){var g=null,i=null;g=dy(a,(1070884>>0));i=ls;document.addEventListener(g,i);}function dy(g,h){var i=null,l=0,j=null;i=String();if((g[h]&255)===0)return String(i);l=0;while(1){j=String.fromCharCode(g[h+l|0]<<24>>24);i=i.concat(j);l=l+1|0;if((g[h+l|0]&255)!==0)continue;break;}return String(i);}function lr(){var n=null,l=null,q=null,p=-0.,o=-0.,i=null,j=null,k=0,g=0;n=fz();l=-64+n|0;dI(l);document.body;if(!(iK|0)){i=dy(a,(1070874>>0));j=document.getElementById(i);iJ=j;iK=1;}i=dy(a,(1070865>>0));q=document.getElementById(i);i=dy(a,(1071000>>0));j=document.getElementById(i);i=dy(a,(1071243>>0));document.getElementById(i);i=dy(a,(1071237>>0));document.getElementById(i);i=q.value;j=j.value;p=+parseInt(i);o=+parseInt(j);f[1049712>>3]=p/1000;f[1049720>>3]=o/1000;hB(hA(hB(1077848|0,1071225|0)|0|0,1049696|0)|0|0,1071221|0)|0;hz();g=32+l|0;hy(g,1049696|0);hA(1077848|0,g)|0;i=iJ;g=16+l|0;lo(g);k=c[8+g>>2];j=a;j=dy(j,k);i.textContent=j;g4(g);g=c[16+(c[1077840>>2]|0)>>2]|0;g=(1077840+g|0);lm(g);g=l|0;ll(g);lk(g);g4(g);dI(n);}function lt(){var i=null,g=0;g=c[16+(c[1077840>>2]|0)>>2]|0;g=(1077840+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1077840>>2]|0)>>2]|0;g=(1077840+g|0);c[8+g>>2]=12;i=fS;+requestAnimationFrame(i);}function fS(){var g=null;lr();g=fS;+requestAnimationFrame(g);}function ls(){var i=null,g=0;g=c[16+(c[1077840>>2]|0)>>2]|0;g=(1077840+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1077840>>2]|0)>>2]|0;g=(1077840+g|0);c[8+g>>2]=12;i=fS;+requestAnimationFrame(i);}function uB(T,U,V){var O=null,n=0,t=0,v=0,w=0,S=0,j=0,x=0,y=0,z=0,A=0,D=0,N=0,P=0,p=-0.,l=0,g=null,G=0,i=-0.,L=0,q=0,o=null;O=fz();g=-173440+O|0;dI(g);n=173360+g|0;hs(n);t=115760+g|0;v=58160+g|0;w=560+g|0;if(!(ih|0)){G=t>>0;o=a;o=o;i=+o.BYTES_PER_ELEMENT;ie[ig]=new Float64Array(o.buffer,(+(uZ(G,~~i)>>>0)),7200);ih=1;}if(!(id|0)){G=v>>0;o=a;o=o;i=+o.BYTES_PER_ELEMENT;ib[ic]=new Float64Array(o.buffer,(+(uZ(G,~~i)>>>0)),7200);id=1;}if(!(ia|0)){G=w>>0;o=a;o=o;i=+o.BYTES_PER_ELEMENT;h_[h$]=new Float64Array(o.buffer,(+(uZ(G,~~i)>>>0)),7200);ia=1;}G=(V|0)===10?1:0;S=V>>>0<2?1:0;j=480+g|0;x=400+g|0;y=320+g|0;z=240+g|0;A=160+g|0;D=80+g|0;N=g|0;L=0;i=-1;while(1){P=uZ(L,200)|0;q=0;while(1){p=(+(q-100|0))/100;f[32+n>>3]=p;f[40+n>>3]=i;f[64+n>>3]=p;f[72+n>>3]=i;hs(j);l=q+P|0;if(S){if((V|0)===0){d0(x,n,U);dZ(j,x);g=t|0;f[(l<<3)+g>>3]=+f[32+j>>3];f[(l<<3)+28800+g>>3]=+f[40+j>>3];f[32+n>>3]=i;f[40+n>>3]=p;f[64+n>>3]=i;f[72+n>>3]=p;d0(A,n,U);dZ(j,A);f[(l<<3)+14400+g>>3]=+f[32+j>>3];f[(l<<3)+43200+g>>3]=+f[40+j>>3];}else{d0(y,n,U);dZ(j,y);g=v|0;f[(l<<3)+g>>3]=+f[64+j>>3];f[(l<<3)+28800+g>>3]=+f[72+j>>3];f[32+n>>3]=i;f[40+n>>3]=p;f[64+n>>3]=i;f[72+n>>3]=p;d0(D,n,U);dZ(j,D);f[(l<<3)+14400+g>>3]=+f[64+j>>3];f[(l<<3)+43200+g>>3]=+f[72+j>>3];}}else if(G){f[64+n>>3]=0;f[72+n>>3]=0;d0(z,n,U);dZ(j,z);g=w|0;f[(l<<3)+g>>3]=+f[32+j>>3];f[(l<<3)+28800+g>>3]=+f[40+j>>3];f[32+n>>3]=i;f[40+n>>3]=p;f[64+n>>3]=0;f[72+n>>3]=0;d0(N,n,U);dZ(j,N);f[(l<<3)+14400+g>>3]=+f[32+j>>3];f[(l<<3)+43200+g>>3]=+f[40+j>>3];}else{f[32+n>>3]=i;f[40+n>>3]=p;f[64+n>>3]=i;f[72+n>>3]=p;}q=q+1|0;if((q|0)!==200)continue;break;}L=L+1|0;if((L|0)===9){g=((V|0)===0?ie:((V|0)===1?ib:h_))[(V|0)===0?ig:((V|0)===1?ic:h$)];dI(O);return g;}i+=.25;continue;}}function uz(n,l,j,i,g){return (l-j)/(i-j);}function uC(t,q,p){var n=null,i=0,o=-0.,j=null,g=0,l=null;n=fz();j=-64+n|0;dI(j);if(!(h9|0)){g=1053104>>0;l=a;l=l;o=+l.BYTES_PER_ELEMENT;h8=new Float64Array(l.buffer,(+(uZ(g,~~o)>>>0)),4);h9=1;}g=32+j|0;k_(g);f[16+g>>3]=q;f[24+g>>3]=p;i=j|0;hy(i,g);f[1053104>>3]=+f[16+i>>3];f[1053112>>3]=+f[24+i>>3];j=h8;dI(n);return j;}function uy(g){}function jC(i,g){i=i|0;g=g|0;k4(a,(i>>0),g);}var iK=0;var iJ=null;var ih=0;var ie=[null];var ig=0;var id=0;var ib=[null];var ic=0;var ia=0;var h_=[null];var h$=0;var h9=0;var h8=null;function uo(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a4;obj.a4.o=2;obj.a4.a=a;return obj;}function un(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a3;obj.a3.o=2;obj.a3.a=a;return obj;}function vj(bytes){var pages=(bytes+65535)>>16;try{__asm.vm.grow(pages);vk(__asm.vm.buffer);return pages<<16;}catch(e){return -1;}}var a=null,c=null,f=null,__asm=null,__heap=null;function vl(){throw new Error('this should be unreachable');};var hB=null;var hA=null;var hz=null;var hy=null;var lo=null;var g4=null;var lm=null;var ll=null;var lk=null;var hs=null;var d0=null;var dZ=null;var k_=null;var fz=null;var dI=null;var Graphics=vl;var JsStruct=vl;vl.promise=vo('derivatives.wasm').then(g=>WebAssembly.instantiate(g,{i:{a3:vl,a2:vl,pT:vl,g5:g5,g8:g8,lt:lt,k5:k5,jC:jC,u$:Math.pow,vj:vj,}})).then(g=>{__asm=g.instance.exports;__heap=__asm.vm.buffer;vk(__heap);hB=__asm.hB;hA=__asm.hA;hz=__asm.hz;hy=__asm.hy;lo=__asm.lo;g4=__asm.g4;lm=__asm.lm;ll=__asm.ll;lk=__asm.lk;hs=__asm.hs;d0=__asm.d0;dZ=__asm.dZ;k_=__asm.k_;fz=__asm.fz;dI=__asm.dI;Graphics=function (){this.i0=0;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}uy(this);};Graphics.prototype.update_array=function (a0,a1){return uC(this,a0,a1);};Graphics.prototype.normalize=function (a0,a1,a2,a3){return uz(this,a0,a1,a2,a3);};Graphics.prototype.conformal_map=function (a0,a1){return uB(this,a0,a1);};Graphics.loadCallback=function (){return ls();};Graphics.initialize=function (){return lt();};Graphics.mainLoop=function (){return lr();};Graphics.init=function (){return uA();};JsStruct=function (a0){this.d0=u0(0.);this.i1=0;this.a2=null;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}uv(this,a0);};JsStruct.prototype.test=function (){return ux(this);};JsStruct.prototype.factorial=function (a0){return uw(this,a0);};Graphics.promise=JsStruct.promise=Promise.resolve();__asm.h6();__asm.h5();__asm.k9();});function vk(g){a=new Uint8Array(g);c=new Int32Array(g);f=new Float64Array(g);}