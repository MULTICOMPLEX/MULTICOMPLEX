"use strict";/*Compiled using Cheerp (R) by Leaning Technologies Ltd*/var xw=Math.imul;var xx=Math.fround;var oSlot=0;var nullArray=[null];var nullObj={d:nullArray,o:0};function xX(p){var b=null,f='function';if(typeof fetch===f)b=fetch(p).then(r=>r.arrayBuffer());else if(typeof require===f){p=require('path').join(__dirname, p);b=new Promise((y,n)=>{require('fs').readFile(p,(e,d)=>{if(e)n(e);else y(d);});});}else b=new Promise((y,n)=>{y(read(p,'binary'));});return b;}function jd(){return +Date.now();}function ja(){return ~~( +Math.random()*4294967296)|0;}function n3(l,m,j){var i=0,g=null;g=n2(l,m,j);i=j-1|0;if((l[m+i|0]&255)===10){g=g.substr(0,i);console.log(g);return;}console.log(g);}function n2(t,u,s){var l=0,n=0,g=null,o=0,p=0,i=0,j=null;g=String();if((s|0)===0)return g;p=s;o=0;while(1){l=t[u+o|0]|0;if((l&255)!==0){n=l&255;if(l<<24>-16777216){i=n;}else if((l&255)<192){i=n&63|i<<6;}else if((l&255)<224){i=n&31;}else if((l&255)<240){i=n&15;}else{i=n&7;}p=p-1|0;o=o+1|0;a:{if((p|0)!==0)if((t[u+o|0]&192)===128)break a;if(i>>>0<65536){j=String.fromCharCode(i);g=g.concat(j);}else{i=i-65536|0;j=String.fromCharCode((i>>>10)+55296|0);g=g.concat(j);j=String.fromCharCode((i&1023)+56320|0);g=g.concat(j);}}if((p|0)!==0)continue;return g;}break;}return g;}function n4(){throw new Error("Abort called");}function w6(l,j){var i=0,g=0;if((j|0)<2)return 1|0;g=1;i=j;while(1){g=xw(g,i)|0;if((i|0)<3)return g|0;i=i-1|0;continue;}}function w7(g){return g.a2;}function w5(n,l){var j=0,i=-0.,g=null;n.d0=l;n.i1=0;j=8393856>>0;g=a;g=g;i=+g.BYTES_PER_ELEMENT;n.a2=new Float64Array(g.buffer,(+(xw(j,~~i)>>>0)),4);}function w_(){var g=null,i=null;g=eG(a,(8411030>>0));i=oB;document.addEventListener(g,i);}function eG(g,h){var i=null,l=0,j=null;i=String();if((g[h]&255)===0)return String(i);l=0;while(1){j=String.fromCharCode(g[h+l|0]<<24>>24);i=i.concat(j);l=l+1|0;if((g[h+l|0]&255)!==0)continue;break;}return String(i);}function oA(){var n=null,l=null,s=null,p=-0.,o=-0.,i=null,j=null,k=0,g=0;n=gP();l=-64+n|0;eY(l);document.body;if(!(k_|0)){i=eG(a,(8411020>>0));j=document.getElementById(i);k9=j;k_=1;}i=eG(a,(8411011>>0));s=document.getElementById(i);i=eG(a,(8411002>>0));j=document.getElementById(i);i=eG(a,(8411144>>0));document.getElementById(i);i=eG(a,(8411462>>0));document.getElementById(i);i=s.value;j=j.value;p=+parseInt(i);o=+parseInt(j);f[8391184>>3]=p/1000;f[8391192>>3]=o/1000;gB(jt(gB(8417272|0,8411450|0)|0|0,8391168|0)|0|0,8411446|0)|0;js();g=32+l|0;hv(g,8391168|0);jt(8417272|0,g)|0;i=k9;g=16+l|0;ox(g);k=c[8+g>>2];j=a;j=eG(j,k);i.textContent=j;i$(g);g=c[16+(c[8417264>>2]|0)>>2]|0;g=(8417264+g|0);ov(g);g=l|0;ou(g);ot(g);i$(g);eY(n);}function oC(){var i=null,g=0;g=c[16+(c[8417264>>2]|0)>>2]|0;g=(8417264+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[8417264>>2]|0)>>2]|0;g=(8417264+g|0);c[8+g>>2]=12;i=g4;+requestAnimationFrame(i);}function g4(){var g=null;oA();g=g4;+requestAnimationFrame(g);}function oB(){var i=null,g=0;g=c[16+(c[8417264>>2]|0)>>2]|0;g=(8417264+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[8417264>>2]|0)>>2]|0;g=(8417264+g|0);c[8+g>>2]=12;i=g4;+requestAnimationFrame(i);}function w$(aP,aS,aV){var aM=null,aE=0,aj=0,af=0,ad=0,ac=0,ab=0,aa=0,$=0,Z=0,al=0,U=-0.,aF=-0.,aG=-0.,aN=0,an=null,aq=null,g=null,h=0,z=0,Q=0,W=0,x=0,p=0,A=0,D=0,s=0,y=0,S=0,R=0,v=0,w=0,O=0,L=0,G=-0.,j=-0.,l=-0.,o=0,T=-0.,V=-0.,au=-0.,ax=-0.,ay=-0.,az=-0.,Y=0,n=0,t=0,i=null;aM=gP();g=-894632+aM|0;eY(g);aE=837032+g|0;aj=779432+g|0;af=721832+g|0;ad=664232+g|0;ac=606632+g|0;ab=549032+g|0;aa=491432+g|0;if(!(kF|0)){o=aE>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;kD[kE]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),7200);kF=1;}$=328216+g|0;Z=165000+g|0;al=1784+g|0;if(!(kC|0)){o=$>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;kA[kB]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),20402);kC=1;}if(!(kz|0)){o=Z>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;kx[ky]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),20402);kz=1;}if(!(kw|0)){o=al>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;ku[kv]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),20402);kw=1;}if(!(kt|0)){o=aj>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;kr[ks]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),7200);kt=1;}if(!(kq|0)){o=af>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;ko[kp]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),7200);kq=1;}if(!(kn|0)){o=ad>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;kl[km]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),7200);kn=1;}if(!(kk|0)){o=ac>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;ki[kj]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),7200);kk=1;}if(!(kh|0)){o=ab>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;kf[kg]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),7200);kh=1;}if(!(ke|0)){o=aa>>0;i=a;i=i;j=+i.BYTES_PER_ELEMENT;kc[kd]=new Float64Array(i.buffer,(+(xw(o,~~j)>>>0)),7200);ke=1;}a:{b:{c:{d:{e:{f:{g:{h:switch(aV|0){case 11:z=1752+g|0;ey(z);Q=1720+g|0;ey(Q);W=1672+g|0;jn(W,2,1);x=1624+g|0;jn(x,3,0);p=1592+g|0;ey(p);A=1560+g|0;ey(A);D=1480+g|0;gw(D);s=1400+g|0;gw(s);y=1392+g|0;oq(y);S=1360+g|0;R=1280+g|0;v=1248+g|0;w=1216+g|0;O=1184+g|0;L=1152+g|0;o=0;l=0;j=0;G=0;while(1){T=0;while(1){f[16+p>>3]=T;f[16+A>>3]=G;fE(D,p);fE(s,A);hm(R,W,D,s);cb(S,R);ex(z,S);oo(v,x,p,A);ex(Q,v);gv(w,z,aS);ex(z,w);gv(O,Q,aS);ex(Q,O);on(L,z);ex(z,L);U=+jl(z);V=+jl(Q);au=+Math.sin(T);au*= +Math.sin(G);ax=+Math.sin(T);ax*= +Math.cos(G);aF=+Math.cos(T);ay=U*au;az=U*ax;aG=U*aF;l+=(aG*aG+(ay*ay+az*az));ay=V*au;az=V*ax;V*=aF;j+=(V*V+(ay*ay+az*az));f[(o<<3)+$>>3]=U*(au*3.1415926535897931);f[(o<<3)+Z>>3]=U*(ax*3.1415926535897931);f[(o<<3)+al>>3]=U*(aF*3.1415926535897931);T+=.031415926535897934;o=o+1|0;if(T<3.1415926535897931)continue;break;}G+=.031415926535897934;if(G<6.2831853071795862)continue;jm(dN(gB(8393304|0,8411437|0)|0|0,l)|0|0);jm(dN(gB(8393304|0,8411428|0)|0|0,j/l)|0|0);break d;}case 10:case 2:case 0:z=(aV|0)!==0?1:0;Q=(aV|0)===10?1:0;W=(aV& -3|0)!==0?1:0;x=1120+g|0;p=1088+g|0;A=1056+g|0;D=1024+g|0;s=944+g|0;y=864+g|0;S=784+g|0;R=704+g|0;v=528+g|0;w=352+g|0;O=176+g|0;L=g|0;o=0;G=-1;i:while(1){j=G*1.5707963267948966;aN=xw(o,200)|0;Y=0;while(1){n=Y+aN|0;l=(+(Y-100|0))/100;if(W){if(Q){ey(x);l*=1.5707963267948966;f[16+x>>3]=l;f[24+x>>3]=j;ey(p);gv(A,x,aS);ex(p,A);g=aE|0;f[(n<<3)+g>>3]=+f[16+p>>3];f[(n<<3)+28800+g>>3]=+f[24+p>>3];f[16+x>>3]=j;f[24+x>>3]=l;gv(D,x,aS);ex(p,D);f[(n<<3)+14400+g>>3]=+f[16+p>>3];f[(n<<3)+43200+g>>3]=+f[24+p>>3];}}else if(z){ji(v);l*=1.5707963267948966;f[48+v>>3]=l;f[56+v>>3]=j;f[80+v>>3]=l;f[88+v>>3]=j;f[128+v>>3]=l;f[136+v>>3]=j;f[160+v>>3]=l;f[168+v>>3]=j;ji(w);jh(O,v,aS);jg(w,O);g=ad|0;f[(n<<3)+g>>3]=+f[48+w>>3];t=n+3600|0;f[(t<<3)+g>>3]=+f[56+w>>3];i=ac|0;f[(n<<3)+i>>3]=+f[80+w>>3];f[(t<<3)+i>>3]=+f[88+w>>3];an=ab|0;f[(n<<3)+an>>3]=+f[128+w>>3];f[(t<<3)+an>>3]=+f[136+w>>3];aq=aa|0;f[(n<<3)+aq>>3]=+f[160+w>>3];f[(t<<3)+aq>>3]=+f[168+w>>3];f[48+v>>3]=j;f[56+v>>3]=l;f[80+v>>3]=j;f[88+v>>3]=l;f[128+v>>3]=j;f[136+v>>3]=l;f[160+v>>3]=j;f[168+v>>3]=l;jh(L,v,aS);jg(w,L);t=n+1800|0;f[(t<<3)+g>>3]=+f[48+w>>3];n=n+5400|0;f[(n<<3)+g>>3]=+f[56+w>>3];f[(t<<3)+i>>3]=+f[80+w>>3];f[(n<<3)+i>>3]=+f[88+w>>3];f[(t<<3)+an>>3]=+f[128+w>>3];f[(n<<3)+an>>3]=+f[136+w>>3];f[(t<<3)+aq>>3]=+f[160+w>>3];f[(n<<3)+aq>>3]=+f[168+w>>3];}else{gw(s);l*=1.5707963267948966;f[32+s>>3]=l;f[40+s>>3]=j;f[64+s>>3]=l;f[72+s>>3]=j;gw(y);jk(S,s,aS);jj(y,S);g=aj|0;f[(n<<3)+g>>3]=+f[32+y>>3];t=n+3600|0;f[(t<<3)+g>>3]=+f[40+y>>3];i=af|0;f[(n<<3)+i>>3]=+f[64+y>>3];f[(t<<3)+i>>3]=+f[72+y>>3];f[32+s>>3]=j;f[40+s>>3]=l;f[64+s>>3]=j;f[72+s>>3]=l;jk(R,s,aS);jj(y,R);t=n+1800|0;f[(t<<3)+g>>3]=+f[32+y>>3];n=n+5400|0;f[(n<<3)+g>>3]=+f[40+y>>3];f[(t<<3)+i>>3]=+f[64+y>>3];f[(n<<3)+i>>3]=+f[72+y>>3];}Y=Y+1|0;if((Y|0)!==200)continue;o=o+1|0;if((o|0)!==9){G+=.25;continue i;}switch(aV|0){case 0:h=ks;g=kr;break a;case 1:break h;case 2:h=km;g=kl;break a;case 3:break g;case 4:break f;case 5:break e;case 10:h=kE;g=kD;break a;case 11:break d;case 12:break c;default:break b;}break;}break;}case 1:break h;case 3:break g;case 4:break f;case 5:break e;case 12:break c;default:break b;}h=kp;g=ko;break a;}h=kj;g=ki;break a;}h=kg;g=kf;break a;}h=kd;g=kc;break a;}h=kB;g=kA;break a;}h=ky;g=kx;break a;}h=kv;g=ku;}g=g[h];eY(aM);return g;}function w9(n,l,j,i,g){return (l-j)/(i-j);}function xa(t,s,p){var n=null,i=0,o=-0.,j=null,g=0,l=null;n=gP();j=-64+n|0;eY(j);if(!(kb|0)){g=8393856>>0;l=a;l=l;o=+l.BYTES_PER_ELEMENT;ka=new Float64Array(l.buffer,(+(xw(g,~~o)>>>0)),4);kb=1;}g=32+j|0;ey(g);f[16+g>>3]=s;f[24+g>>3]=p;i=j|0;hv(i,g);f[8393856>>3]=+f[16+i>>3];f[8393864>>3]=+f[24+i>>3];j=ka;eY(n);return j;}function w8(g){}function lW(i,g){i=i|0;g=g|0;n3(a,(i>>0),g);}var k_=0;var k9=null;var kF=0;var kD=[null];var kE=0;var kC=0;var kA=[null];var kB=0;var kz=0;var kx=[null];var ky=0;var kw=0;var ku=[null];var kv=0;var kt=0;var kr=[null];var ks=0;var kq=0;var ko=[null];var kp=0;var kn=0;var kl=[null];var km=0;var kk=0;var ki=[null];var kj=0;var kh=0;var kf=[null];var kg=0;var ke=0;var kc=[null];var kd=0;var kb=0;var ka=null;function wX(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a4;obj.a4.o=2;obj.a4.a=a;return obj;}function wW(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a3;obj.a3.o=2;obj.a3.a=a;return obj;}function xS(bytes){var pages=(bytes+65535)>>16;try{__asm.xV.grow(pages);xT(__asm.xV.buffer);return pages<<16;}catch(e){return -1;}}var a=null,c=null,f=null,__asm=null,__heap=null;function xU(){throw new Error('this should be unreachable');};var gB=null;var jt=null;var js=null;var hv=null;var ox=null;var i$=null;var ov=null;var ou=null;var ot=null;var ey=null;var jn=null;var gw=null;var oq=null;var dN=null;var jm=null;var fE=null;var hm=null;var cb=null;var ex=null;var oo=null;var gv=null;var on=null;var jl=null;var ar=null;var ai=null;var jk=null;var jj=null;var ji=null;var jh=null;var jg=null;var gP=null;var eY=null;var Graphics=xU;var JsStruct=xU;xU.promise=xX('derivatives.wasm').then(g=>WebAssembly.instantiate(g,{i:{bf:xU,be:xU,nY:xU,ja:ja,jd:jd,oC:oC,n4:n4,lW:lW,xI:Math.pow,xS:xS,}})).then(g=>{__asm=g.instance.exports;__heap=__asm.xV.buffer;xT(__heap);gB=__asm.gB;jt=__asm.jt;js=__asm.js;hv=__asm.hv;ox=__asm.ox;i$=__asm.i$;ov=__asm.ov;ou=__asm.ou;ot=__asm.ot;ey=__asm.ey;jn=__asm.jn;gw=__asm.gw;oq=__asm.oq;dN=__asm.dN;jm=__asm.jm;fE=__asm.fE;hm=__asm.hm;cb=__asm.cb;ex=__asm.ex;oo=__asm.oo;gv=__asm.gv;on=__asm.on;jl=__asm.jl;ar=__asm.ar;ai=__asm.ai;jk=__asm.jk;jj=__asm.jj;ji=__asm.ji;jh=__asm.jh;jg=__asm.jg;gP=__asm.gP;eY=__asm.eY;Graphics=function (){this.i0=0;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}w8(this);};Graphics.prototype.update_array=function (a0,a1){return xa(this,a0,a1);};Graphics.prototype.normalize=function (a0,a1,a2,a3){return w9(this,a0,a1,a2,a3);};Graphics.prototype.conformal_map=function (a0,a1){return w$(this,a0,a1);};Graphics.loadCallback=function (){return oB();};Graphics.initialize=function (){return oC();};Graphics.mainLoop=function (){return oA();};Graphics.init=function (){return w_();};JsStruct=function (a0){this.d0=xx(0.);this.i1=0;this.a2=null;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}w5(this,a0);};JsStruct.prototype.test=function (){return w7(this);};JsStruct.prototype.factorial=function (a0){return w6(this,a0);};Graphics.promise=JsStruct.promise=Promise.resolve();__asm.jz();__asm.jy();__asm.of();});function xT(g){a=new Uint8Array(g);c=new Int32Array(g);f=new Float64Array(g);}