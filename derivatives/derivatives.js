"use strict";/*Compiled using Cheerp (R) by Leaning Technologies Ltd*/var xh=Math.imul;var xi=Math.fround;var oSlot=0;var nullArray=[null];var nullObj={d:nullArray,o:0};function xI(p){var b=null,f='function';if(typeof fetch===f)b=fetch(p).then(r=>r.arrayBuffer());else if(typeof require===f){p=require('path').join(__dirname, p);b=new Promise((y,n)=>{require('fs').readFile(p,(e,d)=>{if(e)n(e);else y(d);});});}else b=new Promise((y,n)=>{y(read(p,'binary'));});return b;}function i8(){return +Date.now();}function i5(){return ~~( +Math.random()*4294967296)|0;}function nU(l,m,j){var i=0,g=null;g=nT(l,m,j);i=j-1|0;if((l[m+i|0]&255)===10){g=g.substr(0,i);console.log(g);return;}console.log(g);}function nT(t,u,q){var l=0,n=0,g=null,o=0,p=0,i=0,j=null;g=String();if((q|0)===0)return g;p=q;o=0;while(1){l=t[u+o|0]|0;if((l&255)!==0){n=l&255;if(l<<24>-16777216){i=n;}else if((l&255)<192){i=n&63|i<<6;}else if((l&255)<224){i=n&31;}else if((l&255)<240){i=n&15;}else{i=n&7;}p=p-1|0;o=o+1|0;a:{if((p|0)!==0)if((t[u+o|0]&192)===128)break a;if(i>>>0<65536){j=String.fromCharCode(i);g=g.concat(j);}else{i=i-65536|0;j=String.fromCharCode((i>>>10)+55296|0);g=g.concat(j);j=String.fromCharCode((i&1023)+56320|0);g=g.concat(j);}}if((p|0)!==0)continue;return g;}break;}return g;}function nh(){throw new Error("Abort called");}function wR(l,j){var i=0,g=0;if((j|0)<2)return 1|0;g=1;i=j;while(1){g=xh(g,i)|0;if((i|0)<3)return g|0;i=i-1|0;continue;}}function wS(g){return g.a2;}function wQ(n,l){var j=0,i=-0.,g=null;n.d0=l;n.i1=0;j=1052904>>0;g=a;g=g;i=+g.BYTES_PER_ELEMENT;n.a2=new Float64Array(g.buffer,(+(xh(j,~~i)>>>0)),4);}function wV(){var g=null,i=null;g=ex(a,(1069922>>0));i=oh;document.addEventListener(g,i);}function ex(g,h){var i=null,l=0,j=null;i=String();if((g[h]&255)===0)return String(i);l=0;while(1){j=String.fromCharCode(g[h+l|0]<<24>>24);i=i.concat(j);l=l+1|0;if((g[h+l|0]&255)!==0)continue;break;}return String(i);}function og(){var n=null,l=null,q=null,p=-0.,o=-0.,i=null,j=null,k=0,g=0;n=gJ();l=-64+n|0;eM(l);document.body;if(!(k2|0)){i=ex(a,(1070548>>0));j=document.getElementById(i);k1=j;k2=1;}i=ex(a,(1070612>>0));q=document.getElementById(i);i=ex(a,(1070693>>0));j=document.getElementById(i);i=ex(a,(1070702>>0));document.getElementById(i);i=ex(a,(1070708>>0));document.getElementById(i);i=q.value;j=j.value;p=+parseInt(i);o=+parseInt(j);f[1050824>>3]=p/1000;f[1050832>>3]=o/1000;jr(jq(jr(1078384|0,1070714|0)|0|0,1050808|0)|0|0,1070745|0)|0;jp();g=32+l|0;ht(g,1050808|0);jq(1078384|0,g)|0;i=k1;g=16+l|0;od(g);k=c[8+g>>2];j=a;j=ex(j,k);i.textContent=j;i4(g);g=c[16+(c[1078376>>2]|0)>>2]|0;g=(1078376+g|0);ob(g);g=l|0;oa(g);n$(g);i4(g);eM(n);}function oi(){var i=null,g=0;g=c[16+(c[1078376>>2]|0)>>2]|0;g=(1078376+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1078376>>2]|0)>>2]|0;g=(1078376+g|0);c[8+g>>2]=12;i=g0;+requestAnimationFrame(i);}function g0(){var g=null;og();g=g0;+requestAnimationFrame(g);}function oh(){var i=null,g=0;g=c[16+(c[1078376>>2]|0)>>2]|0;g=(1078376+g|0);c[4+g>>2]=c[4+g>>2]& -261|4;g=c[16+(c[1078376>>2]|0)>>2]|0;g=(1078376+g|0);c[8+g>>2]=12;i=g0;+requestAnimationFrame(i);}function wW(aK,aN,aR){var aC=null,at=0,aa=0,V=0,U=0,T=0,S=0,W=0,X=0,Z=0,aD=-0.,aA=-0.,aw=-0.,av=-0.,$=0,G=0,D=0,ab=0,ac=0,aF=0,ae=null,af=null,g=null,h=0,t=0,y=0,P=0,Q=0,v=0,w=0,o=-0.,p=0,A=0,n=-0.,l=-0.,aj=-0.,x=0,z=0,L=0,an=0,R=0,j=0,q=0,i=null;aC=gJ();g=-434336+aC|0;eM(g);at=376736+g|0;aa=319136+g|0;V=261536+g|0;U=203936+g|0;T=146336+g|0;S=88736+g|0;W=31136+g|0;if(!(kF|0)){t=at>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;kD[kE]=new Float64Array(i.buffer,(+(xh(t,~~o)>>>0)),7200);kF=1;}if(!(kC|0)){t=aa>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;kA[kB]=new Float64Array(i.buffer,(+(xh(t,~~o)>>>0)),7200);kC=1;}if(!(kz|0)){t=V>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;kx[ky]=new Float64Array(i.buffer,(+(xh(t,~~o)>>>0)),7200);kz=1;}if(!(kw|0)){t=U>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;ku[kv]=new Float64Array(i.buffer,(+(xh(t,~~o)>>>0)),7200);kw=1;}if(!(kt|0)){t=T>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;kr[ks]=new Float64Array(i.buffer,(+(xh(t,~~o)>>>0)),7200);kt=1;}if(!(kq|0)){t=S>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;ko[kp]=new Float64Array(i.buffer,(+(xh(t,~~o)>>>0)),7200);kq=1;}if(!(kn|0)){t=W>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;kl[km]=new Float64Array(i.buffer,(+(xh(t,~~o)>>>0)),7200);kn=1;}t=21200+g|0;X=11264+g|0;Z=1328+g|0;if(!(kk|0)){y=t>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;ki[kj]=new Float64Array(i.buffer,(+(xh(y,~~o)>>>0)),1242);kk=1;}if(!(kh|0)){y=X>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;kf[kg]=new Float64Array(i.buffer,(+(xh(y,~~o)>>>0)),1242);kh=1;}if(!(ke|0)){y=Z>>0;i=a;i=i;o=+i.BYTES_PER_ELEMENT;kc[kd]=new Float64Array(i.buffer,(+(xh(y,~~o)>>>0)),1242);ke=1;}a:{b:{c:{d:{e:{f:{g:{h:switch(aR|0){case 11:y=1296+g|0;fy(y);P=1264+g|0;fy(P);Q=1216+g|0;n9(Q);v=1184+g|0;w=1152+g|0;A=0;p=0;o=-1.5707963267948966;while(1){n=+Math.cos(o);l=+Math.sin(o);aj=+Math.ceil(n*6.2831853071795862/.10134169850289655);x=(A|0)===0||(A|0)===31?1:0;z=~~aj;L=x?1|0:z|0;aD=x?0:6.2831853071795862/(+(z|0));if((L|0)>0){x=L+p|0;l*=9.869604401089358;aj=0;while(1){aA=n* +Math.cos(aj);aw=n* +Math.sin(aj);n7(v,Q,aA,aw);f2(P,v);hn(w,P,aN);f2(y,w);av=+n6(y);f[(p<<3)+t>>3]=aA*9.869604401089358*av;f[(p<<3)+X>>3]=aw*9.869604401089358*av;f[(p<<3)+Z>>3]=l*av;p=p+1|0;if((p|0)!==(x|0)){aj+=aD;continue;}break;}p=x;}A=A+1|0;if((A|0)===32)break d;o+=.10134169850289655;continue;}case 10:case 2:case 0:y=(aR|0)!==0?1:0;P=(aR|0)===10?1:0;Q=(aR& -3|0)!==0?1:0;v=1120+g|0;w=1088+g|0;p=1056+g|0;A=1024+g|0;x=944+g|0;z=864+g|0;L=784+g|0;$=704+g|0;G=528+g|0;D=352+g|0;ab=176+g|0;ac=g|0;an=0;o=-1;i:while(1){n=o*1.5707963267948966;aF=xh(an,200)|0;R=0;while(1){j=R+aF|0;l=(+(R-100|0))/100;if(Q){if(P){fy(v);l*=1.5707963267948966;f[16+v>>3]=l;f[24+v>>3]=n;fy(w);hn(p,v,aN);f2(w,p);g=at|0;f[(j<<3)+g>>3]=+f[16+w>>3];f[(j<<3)+28800+g>>3]=+f[24+w>>3];f[16+v>>3]=n;f[24+v>>3]=l;hn(A,v,aN);f2(w,A);f[(j<<3)+14400+g>>3]=+f[16+w>>3];f[(j<<3)+43200+g>>3]=+f[24+w>>3];}}else if(y){jb(G);l*=1.5707963267948966;f[48+G>>3]=l;f[56+G>>3]=n;f[80+G>>3]=l;f[88+G>>3]=n;f[128+G>>3]=l;f[136+G>>3]=n;f[160+G>>3]=l;f[168+G>>3]=n;jb(D);ja(ab,G,aN);i$(D,ab);g=U|0;f[(j<<3)+g>>3]=+f[48+D>>3];q=j+3600|0;f[(q<<3)+g>>3]=+f[56+D>>3];i=T|0;f[(j<<3)+i>>3]=+f[80+D>>3];f[(q<<3)+i>>3]=+f[88+D>>3];ae=S|0;f[(j<<3)+ae>>3]=+f[128+D>>3];f[(q<<3)+ae>>3]=+f[136+D>>3];af=W|0;f[(j<<3)+af>>3]=+f[160+D>>3];f[(q<<3)+af>>3]=+f[168+D>>3];f[48+G>>3]=n;f[56+G>>3]=l;f[80+G>>3]=n;f[88+G>>3]=l;f[128+G>>3]=n;f[136+G>>3]=l;f[160+G>>3]=n;f[168+G>>3]=l;ja(ac,G,aN);i$(D,ac);q=j+1800|0;f[(q<<3)+g>>3]=+f[48+D>>3];j=j+5400|0;f[(j<<3)+g>>3]=+f[56+D>>3];f[(q<<3)+i>>3]=+f[80+D>>3];f[(j<<3)+i>>3]=+f[88+D>>3];f[(q<<3)+ae>>3]=+f[128+D>>3];f[(j<<3)+ae>>3]=+f[136+D>>3];f[(q<<3)+af>>3]=+f[160+D>>3];f[(j<<3)+af>>3]=+f[168+D>>3];}else{je(x);l*=1.5707963267948966;f[32+x>>3]=l;f[40+x>>3]=n;f[64+x>>3]=l;f[72+x>>3]=n;je(z);jd(L,x,aN);jc(z,L);g=aa|0;f[(j<<3)+g>>3]=+f[32+z>>3];q=j+3600|0;f[(q<<3)+g>>3]=+f[40+z>>3];i=V|0;f[(j<<3)+i>>3]=+f[64+z>>3];f[(q<<3)+i>>3]=+f[72+z>>3];f[32+x>>3]=n;f[40+x>>3]=l;f[64+x>>3]=n;f[72+x>>3]=l;jd($,x,aN);jc(z,$);q=j+1800|0;f[(q<<3)+g>>3]=+f[32+z>>3];j=j+5400|0;f[(j<<3)+g>>3]=+f[40+z>>3];f[(q<<3)+i>>3]=+f[64+z>>3];f[(j<<3)+i>>3]=+f[72+z>>3];}R=R+1|0;if((R|0)!==200)continue;an=an+1|0;if((an|0)!==9){o+=.25;continue i;}switch(aR|0){case 0:h=kB;g=kA;break a;case 1:break h;case 2:h=kv;g=ku;break a;case 3:break g;case 4:break f;case 5:break e;case 10:h=kE;g=kD;break a;case 11:break d;case 12:break c;default:break b;}break;}break;}case 1:break h;case 3:break g;case 4:break f;case 5:break e;case 12:break c;default:break b;}h=ky;g=kx;break a;}h=ks;g=kr;break a;}h=kp;g=ko;break a;}h=km;g=kl;break a;}h=kj;g=ki;break a;}h=kg;g=kf;break a;}h=kd;g=kc;}g=g[h];eM(aC);return g;}function wU(n,l,j,i,g){return (l-j)/(i-j);}function wX(t,q,p){var n=null,i=0,o=-0.,j=null,g=0,l=null;n=gJ();j=-64+n|0;eM(j);if(!(kb|0)){g=1052904>>0;l=a;l=l;o=+l.BYTES_PER_ELEMENT;ka=new Float64Array(l.buffer,(+(xh(g,~~o)>>>0)),4);kb=1;}g=32+j|0;fy(g);f[16+g>>3]=q;f[24+g>>3]=p;i=j|0;ht(i,g);f[1052904>>3]=+f[16+i>>3];f[1052912>>3]=+f[24+i>>3];j=ka;eM(n);return j;}function wT(g){}function lO(i,g){i=i|0;g=g|0;nU(a,(i>>0),g);}var k2=0;var k1=null;var kF=0;var kD=[null];var kE=0;var kC=0;var kA=[null];var kB=0;var kz=0;var kx=[null];var ky=0;var kw=0;var ku=[null];var kv=0;var kt=0;var kr=[null];var ks=0;var kq=0;var ko=[null];var kp=0;var kn=0;var kl=[null];var km=0;var kk=0;var ki=[null];var kj=0;var kh=0;var kf=[null];var kg=0;var ke=0;var kc=[null];var kd=0;var kb=0;var ka=null;function wI(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a4;obj.a4.o=2;obj.a4.a=a;return obj;}function wH(obj){var a=[];a[0]=obj;obj.o=0;obj.a=a;a[1]=obj.a2;obj.a2.o=1;obj.a2.a=a;a[2]=obj.a3;obj.a3.o=2;obj.a3.a=a;return obj;}function xD(bytes){var pages=(bytes+65535)>>16;try{__asm.xG.grow(pages);xE(__asm.xG.buffer);return pages<<16;}catch(e){return -1;}}var a=null,c=null,f=null,__asm=null,__heap=null;function xF(){throw new Error('this should be unreachable');};var jr=null;var jq=null;var jp=null;var ht=null;var od=null;var i4=null;var ob=null;var oa=null;var n$=null;var fy=null;var n9=null;var ak=null;var au=null;var n7=null;var f2=null;var hn=null;var n6=null;var je=null;var jd=null;var jc=null;var jb=null;var ja=null;var i$=null;var gJ=null;var eM=null;var Graphics=xF;var JsStruct=xF;xF.promise=xI('derivatives.wasm').then(g=>WebAssembly.instantiate(g,{i:{ba:xF,a$:xF,nP:xF,i5:i5,i8:i8,oi:oi,nh:nh,lO:lO,xt:Math.pow,xD:xD,}})).then(g=>{__asm=g.instance.exports;__heap=__asm.xG.buffer;xE(__heap);jr=__asm.jr;jq=__asm.jq;jp=__asm.jp;ht=__asm.ht;od=__asm.od;i4=__asm.i4;ob=__asm.ob;oa=__asm.oa;n$=__asm.n$;fy=__asm.fy;n9=__asm.n9;ak=__asm.ak;au=__asm.au;n7=__asm.n7;f2=__asm.f2;hn=__asm.hn;n6=__asm.n6;je=__asm.je;jd=__asm.jd;jc=__asm.jc;jb=__asm.jb;ja=__asm.ja;i$=__asm.i$;gJ=__asm.gJ;eM=__asm.eM;Graphics=function (){this.i0=0;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}wT(this);};Graphics.prototype.update_array=function (a0,a1){return wX(this,a0,a1);};Graphics.prototype.normalize=function (a0,a1,a2,a3){return wU(this,a0,a1,a2,a3);};Graphics.prototype.conformal_map=function (a0,a1){return wW(this,a0,a1);};Graphics.loadCallback=function (){return oh();};Graphics.initialize=function (){return oi();};Graphics.mainLoop=function (){return og();};Graphics.init=function (){return wV();};JsStruct=function (a0){this.d0=xi(0.);this.i1=0;this.a2=null;;this.d=[this];if (arguments.length===1&&arguments[0]===undefined){return;}wQ(this,a0);};JsStruct.prototype.test=function (){return wS(this);};JsStruct.prototype.factorial=function (a0){return wR(this,a0);};Graphics.promise=JsStruct.promise=Promise.resolve();__asm.jx();__asm.jw();__asm.nN();});function xE(g){a=new Uint8Array(g);c=new Int32Array(g);f=new Float64Array(g);}