"use strict";/*Compiled using Cheerp (R) by Leaning Technologies Ltd*/var uq=Math.imul;var ur=Math.fround;var oSlot=0;var nullArray=[null];var nullObj={d:nullArray,o:0};function uR(p){var b=null,f='function';if(typeof fetch===f)b=fetch(p).then(r=>r.arrayBuffer());else if(typeof require===f){p=require('path').join(__dirname, p);b=new Promise((y,n)=>{require('fs').readFile(p,(e,d)=>{if(e)n(e);else y(d);});});}else b=new Promise((y,n)=>{y(read(p,'binary'));});return b;}function f6(){return +Date.now();}function gy(){return ~~( +Math.random()*4294967296)|0;}function lx(l,m,j){var i=0,g=null;g=lw(l,m,j);i=j-1|0;if((l[m+i|0]&255)===10){g=g.substr(0,i);console.log(g);return;}console.log(g);}function lw(r,s,q){var l=0,n=0,g=null,o=0,p=0,i=0,j=null;g=String();if((q|0)===0)return g;p=q;o=0;while(1){l=r[s+o|0]|0;if((l&255)!==0){n=l&255;if(l<<24>-16777216){i=n;}else if((l&255)<192){i=n&63|i<<6;}else if((l&255)<224){i=n&31;}else if((l&255)<240){i=n&15;}else{i=n&7;}p=p-1|0;o=o+1|0;a:{if((p|0)!==0)if((r[s+o|0]&192)===128)break a;if(i>>>0<65536){j=String.fromCharCode(i);g=g.concat(j);}else{i=i-65536|0;j=String.fromCharCode((i>>>10)+55296|0);g=g.concat(j);j=String.fromCharCode((i&1023)+56320|0);g=g.concat(j);}}if((p|0)!==0)continue;return g;}break;}return g;}function ly(){throw new Error("Abort called");}function tZ(l,j){var i=0,g=0;if((j|0)<2)return 1|0;g=1;i=j;while(1){g=uq(g,i)|0;if((i|0)<3)return g|0;i=i-1|0;continue;}}function t0(g){return g.a2;}function tY(n,l){var j=0,i=-0.,g=null;n.d0=l;n.i1=0;j=1052904>>0;g=a;g=g;i=+g.BYTES_PER_ELEMENT;n.a2=new Float64Array(g.buffer,(+(uq(j,~~i)>>>0)),4);}function t3(){var g=null,i=null;g=dJ(a,(1069987>>0));i=k7;document.addEventListener(g,i);}function dJ(g,h){var i=null,l=0,j=null;i=String();if((g[h]&255)===0)return String(i);l=0;while(1){j=String.fromCharCode(g[h+l|0]<<24>>24);i=i.concat(j);l=l+1|0;if((g[h+l|0]&255)!==0)continue;break;}return String(i);}function k6(){var n=null,l=null,q=null,p=-0.,o=-0.,i=null,j=null,k=0,g=0;n=eC();l=-64+n|0;cS(l);document.body;if(!(h8|0)){i=dJ(a,(1070630>>0));j=document.getElementById(i);h7=j;h8=1;}i=dJ(a,(1070569>>0));q=document.getElementById(i);i=dJ(a,(1070738>>0));j=document.getElementById(i);i=dJ(a,(1070747>>0));document.getElementById(i);i=dJ(a,(1070753>>0));document.getElementById(i);i=q.value;j=j.value;p=+parseInt(i);o=+parseInt(j);f[1048856>>3]=p/1000;f[1048864>>3]=o/1000;g7(g6(g7(1077872|0,1070759|0)|0|0,1048840|0)|0|0,1070790|0)|0;g5();g=32+l|0;g4(g,1048840|0);g6(1077872|0,g)|0;i=h7;g=16+l|0;k3(g);k=c[8+g>>2];j=a;j=dJ(j,k);i.textContent=j;gx(g);g=c[16+(c[1077864>>2]|0)>>2]|0;g=(1077864+g|0);k1(g);g=l|0;k0(g);kZ(g);gx(g);cS(n);}function k8(){var i=null,g=0;g=c[16+(c[1077864>>2]|0)>>2]|0;g=(1077864+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1077864>>2]|0)>>2]|0;g=(1077864+g|0);c[8+g>>2]=12;i=fq;+requestAnimationFrame(i);}function fq(){var g=null;k6();g=fq;+requestAnimationFrame(g);}function k7(){var i=null,g=0;g=c[16+(c[1077864>>2]|0)>>2]|0;g=(1077864+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1077864>>2]|0)>>2]|0;g=(1077864+g|0);c[8+g>>2]=12;i=fq;+requestAnimationFrame(i);}function t4(F,I){var z=null,j=0,y=0,n=0,o=null,p=0,q=0,D=0,A=-0.,x=0,i=null,r=null,g=0,v=-0.,w=0,l=0;z=eC();i=-57736+z|0;cS(i);j=57704+i|0;eL(j);y=57696+i|0;kV(y);n=96+i|0;if(!(hM|0)){g=n>>0;r=a;r=r;v=+r.BYTES_PER_ELEMENT;hL=new Float64Array(r.buffer,(+(uq(g,~~v)>>>0)),7200);hM=1;}r=[0];r[0]=I;o=n|0;g=64+i|0;p=32+i|0;q=i|0;w=0;v=-1;while(1){D=uq(w,200)|0;l=0;while(1){A=(+(l-100|0))/100;f[16+j>>3]=A;f[24+j>>3]=v;eL(g);gU(p,r,0,j);gT(g,p);x=l+D|0;f[(x<<3)+o>>3]=+f[16+g>>3];f[(x<<3)+28800+o>>3]=+f[24+g>>3];f[16+j>>3]=v;f[24+j>>3]=A;gU(q,r,0,j);gT(g,q);f[(x<<3)+14400+o>>3]=+f[16+g>>3];f[(x<<3)+43200+o>>3]=+f[24+g>>3];l=l+1|0;if((l|0)!==200)continue;break;}w=w+1|0;if((w|0)===9){i=hL;cS(z);return i;}v+=.25;continue;}}function gU(v,r,s,q){var w=null,j=null,o=0,n=0,g=0,i=0,l=0,p=0;w=eC();j=-6544+w|0;cS(j);o=6512+j|0;eL(o);f[16+o>>3]=0;f[24+o>>3]=1;n=6336+j|0;kR(n);kQ(n,q|0);switch(r[s]|0){case 1:fr(v|0,q|0);break;case 2:g=6304+j|0;fr(g,q|0);em(v|0,g);break;case 3:kP(v|0,q|0);break;case 4:kO(v|0,q|0);break;case 5:kN(v|0,q|0);break;case 6:kM(v|0,q|0);break;case 7:kL(v|0,q|0);break;case 8:kK(v|0,q|0);break;case 9:kJ(v|0,q|0);break;case 10:g=6272+j|0;fr(g,q|0);kI(v|0,g);break;case 11:g=6240+j|0;i=6064+j|0;l=5888+j|0;p=5712+j|0;ej(p,n);dF(l,p);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 12:g=5680+j|0;i=5504+j|0;l=5328+j|0;gS(l,n);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 13:g=5296+j|0;i=5120+j|0;l=4944+j|0;bM(l,n);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 14:g=4912+j|0;i=4736+j|0;l=4560+j|0;bL(l,n);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 15:g=4528+j|0;i=4352+j|0;l=4176+j|0;gR(l,n);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 16:g=4144+j|0;i=3968+j|0;l=3792+j|0;bK(l,n);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 17:g=3760+j|0;i=3584+j|0;l=3408+j|0;bJ(l,n);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 18:g=3376+j|0;i=3200+j|0;l=3024+j|0;gQ(l,n);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 19:g=2992+j|0;i=2816+j|0;l=2640+j|0;p=2464+j|0;ej(p,n);gP(l,p);cp(i,l);bx(g,i);cE(v|0,g,o);break;case 20:g=2288+j|0;i=2112+j|0;ej(i,n);dF(g,i);bx(v|0,g);break;case 21:g=1936+j|0;gS(g,n);bx(v|0,g);break;case 22:g=1760+j|0;bM(g,n);bx(v|0,g);break;case 23:g=1584+j|0;bL(g,n);bx(v|0,g);break;case 24:g=1408+j|0;gR(g,n);bx(v|0,g);break;case 25:g=1232+j|0;bK(g,n);bx(v|0,g);break;case 26:g=1056+j|0;bJ(g,n);bx(v|0,g);break;case 27:g=880+j|0;gQ(g,n);bx(v|0,g);break;case 28:g=704+j|0;i=528+j|0;ej(i,n);gP(g,i);bx(v|0,g);break;default:g=352+j|0;i=176+j|0;l=j|0;ej(l,n);dF(i,l);cp(g,i);bx(v|0,g);}cS(w);}function t2(n,l,j,i,g){return (l-j)/(i-j);}function t5(r,q,p){var n=null,i=0,o=-0.,j=null,g=0,l=null;n=eC();j=-64+n|0;cS(j);if(!(hK|0)){g=1052904>>0;l=a;l=l;o=+l.BYTES_PER_ELEMENT;hJ=new Float64Array(l.buffer,(+(uq(g,~~o)>>>0)),4);hK=1;}g=32+j|0;eL(g);f[16+g>>3]=q;f[24+g>>3]=p;i=j|0;g4(i,g);f[1052904>>3]=+f[16+i>>3];f[1052912>>3]=+f[24+i>>3];j=hJ;cS(n);return j;}function t1(g){}function iD(i,g){i=i|0;g=g|0;lx(a,(i>>0),g);}var h8=0;var h7=null;var hM=0;var hL=null;var hK=0;var hJ=null;function tR(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a4;obj.a4.o=2;obj.a4.a=a;return obj;}function tQ(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a3;obj.a3.o=2;obj.a3.a=a;return obj;}function uM(bytes){var pages=(bytes+65535)>>16;try{__asm.uP.grow(pages);uN(__asm.uP.buffer);return pages<<16;}catch(e){return -1;}}var a=null,c=null,f=null,__asm=null,__heap=null;function uO(){throw new Error('this should be unreachable');};var g7=null;var g6=null;var g5=null;var g4=null;var k3=null;var gx=null;var k1=null;var k0=null;var kZ=null;var eL=null;var kV=null;var gT=null;var kR=null;var kQ=null;var fr=null;var em=null;var kP=null;var kO=null;var kN=null;var kM=null;var kL=null;var kK=null;var kJ=null;var kI=null;var ej=null;var dF=null;var cp=null;var bx=null;var cE=null;var gS=null;var bM=null;var bL=null;var gR=null;var bK=null;var bJ=null;var gQ=null;var gP=null;var eC=null;var cS=null;var Graphics=uO;var JsStruct=uO;uO.promise=uR('derivatives.wasm').then(g=>WebAssembly.instantiate(g,{i:{aJ:uO,aI:uO,lM:uO,gy:gy,f6:f6,k8:k8,ly:ly,iD:iD,uC:Math.pow,uM:uM,}})).then(g=>{__asm=g.instance.exports;__heap=__asm.uP.buffer;uN(__heap);g7=__asm.g7;g6=__asm.g6;g5=__asm.g5;g4=__asm.g4;k3=__asm.k3;gx=__asm.gx;k1=__asm.k1;k0=__asm.k0;kZ=__asm.kZ;eL=__asm.eL;kV=__asm.kV;gT=__asm.gT;kR=__asm.kR;kQ=__asm.kQ;fr=__asm.fr;em=__asm.em;kP=__asm.kP;kO=__asm.kO;kN=__asm.kN;kM=__asm.kM;kL=__asm.kL;kK=__asm.kK;kJ=__asm.kJ;kI=__asm.kI;ej=__asm.ej;dF=__asm.dF;cp=__asm.cp;bx=__asm.bx;cE=__asm.cE;gS=__asm.gS;bM=__asm.bM;bL=__asm.bL;gR=__asm.gR;bK=__asm.bK;bJ=__asm.bJ;gQ=__asm.gQ;gP=__asm.gP;eC=__asm.eC;cS=__asm.cS;Graphics=function (){this.i0=0;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}t1(this);};Graphics.prototype.update_array=function (a0,a1){return t5(this,a0,a1);};Graphics.prototype.normalize=function (a0,a1,a2,a3){return t2(this,a0,a1,a2,a3);};Graphics.prototype.conformal_map=function (a0){return t4(this,a0);};Graphics.loadCallback=function (){return k7();};Graphics.initialize=function (){return k8();};Graphics.mainLoop=function (){return k6();};Graphics.init=function (){return t3();};JsStruct=function (a0){this.d0=ur(0.);this.i1=0;this.a2=null;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}tY(this,a0);};JsStruct.prototype.test=function (){return t0(this);};JsStruct.prototype.factorial=function (a0){return tZ(this,a0);};Graphics.promise=JsStruct.promise=Promise.resolve();__asm.hH();__asm.hG();__asm.lO();});function uN(g){a=new Uint8Array(g);c=new Int32Array(g);f=new Float64Array(g);}