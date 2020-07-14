#include "PM.h"

PM::PM(const Body& body1, const Body& body2, const double& userG): G(userG), Equations(body1, body2)
{
    double o3 = ma*ma;
    double o4=pax*pax;
    double o5=pay*pay;
    double o6=o3+o4+o5;
    double o7=1./sqrt(o6);
    double o16=mb*mb;
    double o17=pbx*pbx;
    double o18=pby*pby;
    double o19=o16+o17+o18;
    double o20=sqrt(o19);
    double o10=o4+o5;
    double o13=1/o6;
    double o21=-qbx;
    double o22=o21+qax;
    double o23=o22*o22;
    double o24=-qby;
    double o25=o24+qay;
    double o26=o25*o25;
    double o27=o23+o26;
    double o28=1./sqrt(o27);
    double o9=sqrt(o6);
    double o11=pow(o6,-2.);
    double o12=-2.*o10*o11*pax;
    double o14=2.*o13*pax;
    double o15=o12+o14;
    double o30=o10*o13;
    double o31=o17+o18;
    double o32=1/o19;
    double o33=o31*o32;
    double o34=1.+o30+o33;
    double o36=-qax;
    double o37=o36+qbx;
    double o38=o37*o37;
    double o39=-qay;
    double o40=o39+qby;
    double o41=o40*o40;
    double o42=o38+o41;
    double o43=1./sqrt(o42);
    double o48=7.*pbx;
    double o73=pax*pbx;
    double o74=pay*pby;
    double o75=o73+o74;
    double o76=o75*o75;
    double o64=o22*o28*pax;
    double o65=o25*o28*pay;
    double o66=o64+o65;
    double o67=o66*o66;
    double o49=o22*o28*pbx;
    double o50=o25*o28*pby;
    double o51=o49+o50;
    double o68=o3+o67;
    double o69=sqrt(o68);
    double o81=o51*o51;
    double o89=o67*o81;
    double o101=1./sqrt(o68);
    double o63=1./sqrt(o19);
    double o70=o69*o7;
    double o71=1.+o70;
    double o77=-(o10*o76);
    double o78=2.*o67*o76;
    double o79=-2.*o10*o51*o66*o75;
    double o80=o10*o10;
    double o82=o80*o81;
    double o83=o77+o78+o79+o82;
    double o102=pow(o6,-1.5);
    double o85=o10*o31;
    double o86=-3.*o31*o67;
    double o87=8.*o51*o66*o75;
    double o88=-3.*o10*o81;
    double o90=o85+o86+o87+o88+o89;
    double o92=-o17;
    double o93=-o18;
    double o94=o92+o93;
    double o125=2.*o22*o28*o66*o81;
    double o107=pow(o71,-2.);
    double o84=2.*o13*o83;
    double o91=o69*o7*o90;
    double o95=o67*o94;
    double o96=2.*o51*o66*o75;
    double o97=-(o10*o81);
    double o98=o76+o89+o95+o96+o97;
    double o99=2.*o98;
    double o100=o84+o91+o99;
    double o55=o37*o43*pbx;
    double o56=o40*o43*pby;
    double o57=o55+o56;
    double o140=o57*o57;
    double o141=o140+o16;
    double o149=o37*o43*pax;
    double o150=o40*o43*pay;
    double o151=o149+o150;
    double o143=sqrt(o141);
    double o128=2.*o31*pax;
    double o120=2.*o75*pbx;
    double o162=2.*o140*o151*o37*o43;
    double o142=1./sqrt(o141);
    double o144=o143*o63;
    double o145=1.+o144;
    double o146=pow(o145,-2.);
    double o148=o31*o31;
    double o174=o151*o151;
    double o186=o140*o174;
    double o200=-2.*o10*o11*pay;
    double o201=2.*o13*pay;
    double o202=o200+o201;
    double o209=7.*pby;
    double o72=pow(o71,-3.);
    double o239=2.*o25*o28*o66*o81;
    double o138=pow(o68,-1.5);
    double o242=2.*o31*pay;
    double o234=2.*o75*pby;
    double o264=2.*o140*o151*o40*o43;
    double o173=-(o31*o76);
    double o175=o148*o174;
    double o176=-2.*o151*o31*o57*o75;
    double o177=2.*o140*o76;
    double o178=o173+o175+o176+o177;
    double o179=2.*o178*o32;
    double o180=-(o174*o31);
    double o181=2.*o151*o57*o75;
    double o182=-o4;
    double o183=-o5;
    double o184=o182+o183;
    double o185=o140*o184;
    double o187=o180+o181+o185+o186+o76;
    double o188=2.*o187;
    double o189=-3.*o174*o31;
    double o190=8.*o151*o57*o75;
    double o191=-3.*o10*o140;
    double o192=o186+o189+o190+o191+o85;
    double o193=o143*o192*o63;
    double o194=o179+o188+o193;
    double o280=pow(o19,-2.);
    double o281=-2.*o280*o31*pbx;
    double o282=2.*o32*pbx;
    double o283=o281+o282;
    double o290=7.*pax;
    double o313=2.*o22*o28*o51*o67;
    double o316=2.*o75*pax;
    double o308=2.*o10*pbx;
    double o329=2.*o174*o37*o43*o57;
    double o299=pow(o19,-1.5);
    double o364=-2.*o280*o31*pby;
    double o365=2.*o32*pby;
    double o366=o364+o365;
    double o373=7.*pay;
    double o395=2.*o25*o28*o51*o67;
    double o398=2.*o75*pay;
    double o390=2.*o10*pby;
    double o411=2.*o174*o40*o43*o57;
    double o356=pow(o145,-3.);
    double o358=pow(o141,-1.5);
    
    qdot0[0][0]=0.25*G*(o28*(o48+o22*o28*o51)+o43*(o48+o37*o43*o57))+o7*pax-0.5*G*(o15*o20*o28*o9+o15*o20*o43*o9+o20*o28*o34*o7*pax+o20*o34*o43*o7*pax)-0.25*G*(-(o100*o107*o138*o22*o28*o43*o63*o66)-o102*o142*o146*o194*o28*pax-2.*o100*o101*o43*o63*o72*(o101*o22*o28*o66*o7-o102*o69*pax)+o142*o146*o28*o7*(2.*(o120+o162-2.*o151*o31*o37*o43+2.*o37*o43*o57*o75-2.*o140*pax+2.*o151*o57*pbx)+o143*o63*(o128+o162-6.*o151*o31*o37*o43+8.*o37*o43*o57*o75-6.*o140*pax+8.*o151*o57*pbx)+2.*o32*(2.*o148*o151*o37*o43-2.*o31*o37*o43*o57*o75-2.*o151*o31*o57*pbx+4.*o140*o75*pbx-2.*o31*o75*pbx))+o101*o107*o43*o63*(o101*o22*o28*o66*o7*o90-4.*o11*o83*pax-o102*o69*o90*pax+2.*(o120+o125+2.*o22*o28*o51*o75+2.*o22*o28*o66*o94-2.*o81*pax+2.*o51*o66*pbx)+o69*o7*(o125+o128-6.*o22*o28*o31*o66+8.*o22*o28*o51*o75-6.*o81*pax+8.*o51*o66*pbx)+2.*o13*(-2.*o10*o22*o28*o51*o75+4.*o22*o28*o66*o76-4.*o51*o66*o75*pax-2.*o76*pax+4.*o10*o81*pax-2.*o10*o51*o66*pbx-2.*o10*o75*pbx+4.*o67*o75*pbx)));
    
    qdot0[0][1]=0.25*G*(o28*(o209+o25*o28*o51)+o43*(o209+o40*o43*o57))+o7*pay-0.5*G*(o20*o202*o28*o9+o20*o202*o43*o9+o20*o28*o34*o7*pay+o20*o34*o43*o7*pay)-0.25*G*(-(o100*o107*o138*o25*o28*o43*o63*o66)-o102*o142*o146*o194*o28*pay-2.*o100*o101*o43*o63*o72*(o101*o25*o28*o66*o7-o102*o69*pay)+o142*o146*o28*o7*(2.*(o234+o264-2.*o151*o31*o40*o43+2.*o40*o43*o57*o75-2.*o140*pay+2.*o151*o57*pby)+o143*o63*(o242+o264-6.*o151*o31*o40*o43+8.*o40*o43*o57*o75-6.*o140*pay+8.*o151*o57*pby)+2.*o32*(2.*o148*o151*o40*o43-2.*o31*o40*o43*o57*o75-2.*o151*o31*o57*pby+4.*o140*o75*pby-2.*o31*o75*pby))+o101*o107*o43*o63*(o101*o25*o28*o66*o7*o90-4.*o11*o83*pay-o102*o69*o90*pay+2.*(o234+o239+2.*o25*o28*o51*o75+2.*o25*o28*o66*o94-2.*o81*pay+2.*o51*o66*pby)+o69*o7*(o239+o242-6.*o25*o28*o31*o66+8.*o25*o28*o51*o75-6.*o81*pay+8.*o51*o66*pby)+2.*o13*(-2.*o10*o25*o28*o51*o75+4.*o25*o28*o66*o76-4.*o51*o66*o75*pay-2.*o76*pay+4.*o10*o81*pay-2.*o10*o51*o66*pby-2.*o10*o75*pby+4.*o67*o75*pby)));
    
    qdot0[1][0]=0.25*G*(o43*(o290+o151*o37*o43)+o28*(o290+o22*o28*o66))+o63*pbx-0.5*G*(o20*o28*o283*o9+o20*o283*o43*o9+o28*o34*o63*o9*pbx+o34*o43*o63*o9*pbx)-0.25*G*(-(o146*o194*o28*o358*o37*o43*o57*o7)-o100*o101*o107*o299*o43*pbx-2.*o142*o194*o28*o356*o7*(o142*o37*o43*o57*o63-o143*o299*pbx)+o101*o107*o43*o63*(2.*o13*(-2.*o10*o22*o28*o66*o75+2.*o22*o28*o51*o80-2.*o10*o51*o66*pax-2.*o10*o75*pax+4.*o67*o75*pax)+o69*o7*(o308+o313-6.*o10*o22*o28*o51+8.*o22*o28*o66*o75+8.*o51*o66*pax-6.*o67*pbx)+2.*(o313+o316-2.*o10*o22*o28*o51+2.*o22*o28*o66*o75+2.*o51*o66*pax-2.*o67*pbx))+o142*o146*o28*o7*(o142*o192*o37*o43*o57*o63-4.*o178*o280*pbx-o143*o192*o299*pbx+o143*o63*(o308+o329-6.*o10*o37*o43*o57+8.*o151*o37*o43*o75+8.*o151*o57*pax-6.*o174*pbx)+2.*(o316+o329+2.*o184*o37*o43*o57+2.*o151*o37*o43*o75+2.*o151*o57*pax-2.*o174*pbx)+2.*o32*(-2.*o151*o31*o37*o43*o75+4.*o37*o43*o57*o76-2.*o151*o31*o57*pax+4.*o140*o75*pax-2.*o31*o75*pax+4.*o174*o31*pbx-4.*o151*o57*o75*pbx-2.*o76*pbx)));
    
    qdot0[1][1]=0.25*G*(o43*(o373+o151*o40*o43)+o28*(o373+o25*o28*o66))+o63*pby-0.5*G*(o20*o28*o366*o9+o20*o366*o43*o9+o28*o34*o63*o9*pby+o34*o43*o63*o9*pby)-0.25*G*(-(o146*o194*o28*o358*o40*o43*o57*o7)-o100*o101*o107*o299*o43*pby-2.*o142*o194*o28*o356*o7*(o142*o40*o43*o57*o63-o143*o299*pby)+o101*o107*o43*o63*(2.*o13*(-2.*o10*o25*o28*o66*o75+2.*o25*o28*o51*o80-2.*o10*o51*o66*pay-2.*o10*o75*pay+4.*o67*o75*pay)+o69*o7*(o390+o395-6.*o10*o25*o28*o51+8.*o25*o28*o66*o75+8.*o51*o66*pay-6.*o67*pby)+2.*(o395+o398-2.*o10*o25*o28*o51+2.*o25*o28*o66*o75+2.*o51*o66*pay-2.*o67*pby))+o142*o146*o28*o7*(o142*o192*o40*o43*o57*o63-4.*o178*o280*pby-o143*o192*o299*pby+o143*o63*(o390+o411-6.*o10*o40*o43*o57+8.*o151*o40*o43*o75+8.*o151*o57*pay-6.*o174*pby)+2.*(o398+o411+2.*o184*o40*o43*o57+2.*o151*o40*o43*o75+2.*o151*o57*pay-2.*o174*pby)+2.*o32*(-2.*o151*o31*o40*o43*o75+4.*o40*o43*o57*o76-2.*o151*o31*o57*pay+4.*o140*o75*pay-2.*o31*o75*pay+4.*o174*o31*pby-4.*o151*o57*o75*pby-2.*o76*pby)));
    
    /***********************/
    
    double z6=pax*pax;
    double z7=pay*pay;
    double z5=ma*ma;
    double z8=z5+z6+z7;
    double z11=pbx*pbx;
    double z13=pby*pby;
    double z10=mb*mb;
    double z14=z10+z11+z13;
    double z23=-qbx;
    double z24=qax+z23;
    double z9=sqrt(z8);
    double z15=sqrt(z14);
    double z16=z6+z7;
    double z17=1/z8;
    double z18=z16*z17;
    double z19=z11+z13;
    double z20=1/z14;
    double z21=z19*z20;
    double z22=1.+z18+z21;
    double z33=-qax;
    double z34=qbx+z33;
    double z25=z24*z24;
    double z26=-qby;
    double z27=qay+z26;
    double z28=z27*z27;
    double z29=z25+z28;
    double z48=1./sqrt(z29);
    double z30=pow(z29,-1.5);
    double z49=pax*z24*z48;
    double z50=pay*z27*z48;
    double z51=z49+z50;
    double z52=pbx*z24*z48;
    double z53=pby*z27*z48;
    double z54=z52+z53;
    double z35=z34*z34;
    double z36=-qay;
    double z37=qby+z36;
    double z38=z37*z37;
    double z39=z35+z38;
    double z40=pow(z39,-1.5);
    double z70=1./sqrt(z39);
    double z44=pax*pbx;
    double z45=pay*pby;
    double z46=z44+z45;
    double z47=7.*z46;
    double z75=pax*z34*z70;
    double z76=pay*z37*z70;
    double z77=z75+z76;
    double z83=pbx*z34*z70;
    double z84=pby*z37*z70;
    double z85=z83+z84;
    double z102=z46*z46;
    double z96=z51*z51;
    double z95=1./sqrt(z8);
    double z97=z5+z96;
    double z98=sqrt(z97);
    double z107=z54*z54;
    double z115=z107*z96;
    double z94=1./sqrt(z14);
    double z99=z95*z98;
    double z100=1.+z99;
    double z101=pow(z100,-2.);
    double z63=-(pax*z25*z30);
    double z64=pax*z48;
    double z65=-(pay*z24*z27*z30);
    double z66=z63+z64+z65;
    double z58=-(pbx*z25*z30);
    double z59=pbx*z48;
    double z60=-(pby*z24*z27*z30);
    double z61=z58+z59+z60;
    double z106=z16*z16;
    double z118=-z11;
    double z119=-z13;
    double z120=z118+z119;
    double z139=2.*z54*z61*z96;
    double z140=2.*z107*z51*z66;
    double z127=1./sqrt(z97);
    double z111=z16*z19;
    double z112=-3.*z19*z96;
    double z113=8.*z46*z51*z54;
    double z114=-3.*z107*z16;
    double z116=z111+z112+z113+z114+z115;
    double z103=-(z102*z16);
    double z104=2.*z102*z96;
    double z105=-2.*z16*z46*z51*z54;
    double z108=z106*z107;
    double z109=z103+z104+z105+z108;
    double z110=2.*z109*z17;
    double z117=z116*z95*z98;
    double z121=z120*z96;
    double z122=2.*z46*z51*z54;
    double z123=-(z107*z16);
    double z124=z102+z115+z121+z122+z123;
    double z125=2.*z124;
    double z126=z110+z117+z125;
    double z157=z85*z85;
    double z158=z10+z157;
    double z79=pax*z35*z40;
    double z80=pay*z34*z37*z40;
    double z81=-(pax*z70);
    double z82=z79+z80+z81;
    double z71=pbx*z35*z40;
    double z72=pby*z34*z37*z40;
    double z73=-(pbx*z70);
    double z74=z71+z72+z73;
    double z160=sqrt(z158);
    double z178=z77*z77;
    double z179=2.*z178*z74*z85;
    double z180=2.*z157*z77*z82;
    double z159=1./sqrt(z158);
    double z161=z160*z94;
    double z162=1.+z161;
    double z164=z19*z19;
    double z174=-z6;
    double z175=-z7;
    double z176=z174+z175;
    double z192=z157*z178;
    double z189=-3.*z178*z19;
    double z190=8.*z46*z77*z85;
    double z191=-3.*z157*z16;
    double z193=z111+z189+z190+z191+z192;
    double z163=pow(z162,-2.);
    double z199=-(z102*z19);
    double z200=z164*z178;
    double z201=-2.*z19*z46*z77*z85;
    double z202=2.*z102*z157;
    double z203=z199+z200+z201+z202;
    double z204=2.*z20*z203;
    double z205=-(z178*z19);
    double z206=2.*z46*z77*z85;
    double z207=z157*z176;
    double z208=z102+z192+z205+z206+z207;
    double z209=2.*z208;
    double z210=z160*z193*z94;
    double z211=z204+z209+z210;
    double z55=z51*z54;
    double z56=z47+z55;
    double z89=z77*z85;
    double z90=z47+z89;
    double z223=pay*z48;
    double z224=-(pax*z24*z27*z30);
    double z225=-(pay*z28*z30);
    double z226=z223+z224+z225;
    double z228=pby*z48;
    double z229=-(pbx*z24*z27*z30);
    double z230=-(pby*z28*z30);
    double z231=z228+z229+z230;
    double z260=2.*z107*z226*z51;
    double z263=2.*z231*z54*z96;
    double z152=pow(z97,-1.5);
    double z154=pow(z100,-3.);
    double z155=1/z97;
    double z241=pax*z34*z37*z40;
    double z242=pay*z38*z40;
    double z243=-(pay*z70);
    double z244=z241+z242+z243;
    double z236=pbx*z34*z37*z40;
    double z237=pby*z38*z40;
    double z238=-(pby*z70);
    double z239=z236+z237+z238;
    double z287=2.*z178*z239*z85;
    double z288=2.*z157*z244*z77;
    double z197=1/z158;
    double z198=pow(z162,-3.);
    double z213=pow(z158,-1.5);
    double z316=pax*z25*z30;
    double z317=-(pax*z48);
    double z318=pay*z24*z27*z30;
    double z319=z316+z317+z318;
    double z311=pbx*z25*z30;
    double z312=-(pbx*z48);
    double z313=pby*z24*z27*z30;
    double z314=z311+z312+z313;
    double z349=2.*z314*z54*z96;
    double z350=2.*z107*z319*z51;
    double z328=-(pax*z35*z40);
    double z329=-(pay*z34*z37*z40);
    double z330=pax*z70;
    double z331=z328+z329+z330;
    double z323=-(pbx*z35*z40);
    double z324=-(pby*z34*z37*z40);
    double z325=pbx*z70;
    double z326=z323+z324+z325;
    double z374=2.*z178*z326*z85;
    double z375=2.*z157*z331*z77;
    double z397=-(pay*z48);
    double z398=pax*z24*z27*z30;
    double z399=pay*z28*z30;
    double z400=z397+z398+z399;
    double z402=-(pby*z48);
    double z403=pbx*z24*z27*z30;
    double z404=pby*z28*z30;
    double z405=z402+z403+z404;
    double z434=2.*z107*z400*z51;
    double z437=2.*z405*z54*z96;
    double z415=-(pax*z34*z37*z40);
    double z416=-(pay*z38*z40);
    double z417=pay*z70;
    double z418=z415+z416+z417;
    double z410=-(pbx*z34*z37*z40);
    double z411=-(pby*z38*z40);
    double z412=pby*z70;
    double z413=z410+z411+z412;
    double z461=2.*z178*z413*z85;
    double z462=2.*z157*z418*z77;
    
    pdot0[0][0]=0.5*G*(-(z15*z22*z24*z30*z9)+z15*z22*z34*z40*z9)-0.25*G*(-(z24*z30*z56)+z48*(z51*z61+z54*z66)+z70*(z74*z77+z82*z85)+z34*z40*z90)+0.25*G*(z101*z126*z127*z34*z40*z94-z101*z126*z152*z51*z66*z70*z94-z159*z163*z211*z24*z30*z95-z163*z211*z213*z48*z74*z85*z95-2.*z126*z154*z155*z51*z66*z70*z94*z95-2.*z197*z198*z211*z48*z74*z85*z94*z95+z159*z163*z48*(2.*(z179+z180+2.*z46*z74*z77-2.*z19*z77*z82+2.*z176*z74*z85+2.*z46*z82*z85)+2.*z20*(-2.*z19*z46*z74*z77+2.*z164*z77*z82+4.*z102*z74*z85-2.*z19*z46*z82*z85)+z159*z193*z74*z85*z94+z160*(z179+z180+8.*z46*z74*z77-6.*z19*z77*z82-6.*z16*z74*z85+8.*z46*z82*z85)*z94)*z95+z101*z127*z70*z94*(2.*(z139+z140+2.*z46*z51*z61-2.*z16*z54*z61+2.*z120*z51*z66+2.*z46*z54*z66)+2.*z17*(-2.*z16*z46*z51*z61+2.*z106*z54*z61+4.*z102*z51*z66-2.*z16*z46*z54*z66)+z116*z127*z51*z66*z95+(z139+z140+8.*z46*z51*z61-6.*z16*z54*z61-6.*z19*z51*z66+8.*z46*z54*z66)*z95*z98));
    
    pdot0[0][1]=0.5*G*(-(z15*z22*z27*z30*z9)+z15*z22*z37*z40*z9)-0.25*G*(z48*(z231*z51+z226*z54)-z27*z30*z56+z70*(z239*z77+z244*z85)+z37*z40*z90)+0.25*G*(z101*z126*z127*z37*z40*z94-z101*z126*z152*z226*z51*z70*z94-z159*z163*z211*z27*z30*z95-z163*z211*z213*z239*z48*z85*z95-2.*z126*z154*z155*z226*z51*z70*z94*z95-2.*z197*z198*z211*z239*z48*z85*z94*z95+z159*z163*z48*(2.*(z287+z288-2.*z19*z244*z77+2.*z239*z46*z77+2.*z176*z239*z85+2.*z244*z46*z85)+2.*z20*(2.*z164*z244*z77-2.*z19*z239*z46*z77+4.*z102*z239*z85-2.*z19*z244*z46*z85)+z159*z193*z239*z85*z94+z160*(z287+z288-6.*z19*z244*z77+8.*z239*z46*z77-6.*z16*z239*z85+8.*z244*z46*z85)*z94)*z95+z101*z127*z70*z94*(2.*(z260+z263+2.*z120*z226*z51+2.*z231*z46*z51-2.*z16*z231*z54+2.*z226*z46*z54)+2.*z17*(4.*z102*z226*z51-2.*z16*z231*z46*z51+2.*z106*z231*z54-2.*z16*z226*z46*z54)+z116*z127*z226*z51*z95+(z260+z263-6.*z19*z226*z51+8.*z231*z46*z51-6.*z16*z231*z54+8.*z226*z46*z54)*z95*z98));
    
    pdot0[1][0]=0.5*G*(z15*z22*z24*z30*z9-z15*z22*z34*z40*z9)-0.25*G*(z48*(z314*z51+z319*z54)+z24*z30*z56+z70*(z326*z77+z331*z85)-z34*z40*z90)+0.25*G*(-(z101*z126*z127*z34*z40*z94)-z101*z126*z152*z319*z51*z70*z94+z159*z163*z211*z24*z30*z95-z163*z211*z213*z326*z48*z85*z95-2.*z126*z154*z155*z319*z51*z70*z94*z95-2.*z197*z198*z211*z326*z48*z85*z94*z95+z159*z163*z48*(2.*(z374+z375-2.*z19*z331*z77+2.*z326*z46*z77+2.*z176*z326*z85+2.*z331*z46*z85)+2.*z20*(2.*z164*z331*z77-2.*z19*z326*z46*z77+4.*z102*z326*z85-2.*z19*z331*z46*z85)+z159*z193*z326*z85*z94+z160*(z374+z375-6.*z19*z331*z77+8.*z326*z46*z77-6.*z16*z326*z85+8.*z331*z46*z85)*z94)*z95+z101*z127*z70*z94*(2.*(z349+z350+2.*z120*z319*z51+2.*z314*z46*z51-2.*z16*z314*z54+2.*z319*z46*z54)+2.*z17*(4.*z102*z319*z51-2.*z16*z314*z46*z51+2.*z106*z314*z54-2.*z16*z319*z46*z54)+z116*z127*z319*z51*z95+(z349+z350-6.*z19*z319*z51+8.*z314*z46*z51-6.*z16*z314*z54+8.*z319*z46*z54)*z95*z98));
    
    pdot0[1][1]=0.5*G*(z15*z22*z27*z30*z9-z15*z22*z37*z40*z9)-0.25*G*(z48*(z405*z51+z400*z54)+z27*z30*z56+z70*(z413*z77+z418*z85)-z37*z40*z90)+0.25*G*(-(z101*z126*z127*z37*z40*z94)-z101*z126*z152*z400*z51*z70*z94+z159*z163*z211*z27*z30*z95-z163*z211*z213*z413*z48*z85*z95-2.*z126*z154*z155*z400*z51*z70*z94*z95-2.*z197*z198*z211*z413*z48*z85*z94*z95+z159*z163*z48*(2.*(z461+z462-2.*z19*z418*z77+2.*z413*z46*z77+2.*z176*z413*z85+2.*z418*z46*z85)+2.*z20*(2.*z164*z418*z77-2.*z19*z413*z46*z77+4.*z102*z413*z85-2.*z19*z418*z46*z85)+z159*z193*z413*z85*z94+z160*(z461+z462-6.*z19*z418*z77+8.*z413*z46*z77-6.*z16*z413*z85+8.*z418*z46*z85)*z94)*z95+z101*z127*z70*z94*(2.*(z434+z437+2.*z120*z400*z51+2.*z405*z46*z51-2.*z16*z405*z54+2.*z400*z46*z54)+2.*z17*(4.*z102*z400*z51-2.*z16*z405*z46*z51+2.*z106*z405*z54-2.*z16*z400*z46*z54)+z116*z127*z400*z51*z95+(z434+z437-6.*z19*z400*z51+8.*z405*z46*z51-6.*z16*z405*z54+8.*z400*z46*z54)*z95*z98));
}



