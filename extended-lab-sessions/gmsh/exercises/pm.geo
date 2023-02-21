Include "pm_data.pro";

Function Circ_
  pcirc()={}; lcirc()={};
  pcirc()+=newp; Point(pcirc(0)) = {x0+radius, y0,        0., lc_circ};
  pcirc()+=newp; Point(pcirc(1)) = {x0,        y0+radius, 0., lc_circ};
  pcirc()+=newp; Point(pcirc(2)) = {x0-radius, y0,        0., lc_circ};
  pcirc()+=newp; Point(pcirc(3)) = {x0,        y0-radius, 0., lc_circ};

  For k In {0:3}
    k_ = (k<3)?k+1:0;
    lcirc()+=newl; Circle(lcirc(k))  = {pcirc(k),ccirc,pcirc(k_)};
  EndFor

  // Line Loop(newll) = {lcirc()};
  // Plane Surface(news) = {newll-1};
Return




x0=0.; y0=0.;
cen0 = newp; Point(newp)={x0, y0, 0, lcs*2};

pStatorSlotX1() = {p1x, p2x, p3x, p4x, p5x};
pStatorSlotY1() = {p1y, p2y, p3y, p4y, p5y};

pStatorSlotX2() = {p5x, p4x, p8x, p7x, p6x};
pStatorSlotY2() = {p5y, p4y, p8y, p7y, p6y};

pStatorSlotX3() = {p1x, p5x, p6x, p10x, p9x};
pStatorSlotY3() = {p1y, p5y, p6y, p10y, p9y};

For k In {0:4}
  pnts1()+=newp; Point(newp)={pStatorSlotX1(k), pStatorSlotY1(k) , 0, lcs};
  pnts2()+=newp; Point(newp)={pStatorSlotX2(k), pStatorSlotY2(k) , 0, lcs};
  pnts3()+=newp; Point(newp)={pStatorSlotX3(k), pStatorSlotY3(k) , 0, lcs};
EndFor

For k In {0:#pnts1()-1}
  k2 = (k<#pnts1()-1)?k+1:0;
  lin1()+=newl; Line(newl) = {pnts1(k), pnts1(k2)};
  lin2()+=newl; Line(newl) = {pnts2(k), pnts2(k2)};
  lin3()+=newl; Line(newl) = {pnts3(k), pnts3(k2)};
EndFor

Line Loop(newll) = {lin1()};
surf_slot()+=news; Plane Surface(news) = {newll-1};
Line Loop(newll) = {lin2()};
surf_slot()+=news; Plane Surface(news) = {newll-1};
Line Loop(newll) = {lin3()};
surf_slot()+=news; Plane Surface(news) = {newll-1};


For ind In {1:(Ns-1)}
  surf_slot() += Rotate {{0, 0, 1}, {0, 0, 0}, 2*ind*Pi/Ns} { Duplicata{Surface{surf_slot({0:2})};} };
EndFor


// pRotorOpenX() = {p11x, p12x, p13x, p14x};
// pRotorOpenY() = {p11y, p12y, p13y, p14y};

// For k In {0:3}
//   pntr()+=newp; Point(newp)={pRotorOpenX(k), pRotorOpenY(k), 0, lcs};
// EndFor


ccirc = cen0;

radius =Dshaft/2;
lc_circ = lcs*2;
Call Circ_ ;
llshaft=newll; Line Loop(newll) = {lcirc()};
surf_shaft=news; Plane Surface(surf_shaft) = {newll-1};


radius =  D2s/2-lag-Hpm;
lc_circ = lcs*2;
Call Circ_ ;
llring=newll; Line Loop(newll) = {lcirc()};
surf_ring=news; Plane Surface(surf_ring) = {newll-1, llshaft};

radius =  D2s/2-lag;
lc_circ = lcs;
Call Circ_ ;
llring2=newll; Line Loop(newll) = {lcirc()};
surf_ring2=news; Plane Surface(surf_ring2) = {newll-1, llring};

pnts_close_slot0() = {13, 66, 120, 174, 228, 282, 336, 390, 444, 498, 552, 606};
pnts_close_slot1() = {16, 70, 124, 178, 232, 286, 340, 394, 448, 502, 556, 610};

For k In {0:Ns-1}
  k_ = (k<Ns-1)?k+1:0;
  Circle(newl) = {pnts_close_slot1(k), cen0, pnts_close_slot0(k_)};
EndFor

llslots=newll;
Curve Loop(newll) = {248, 218, 249, 12, 238, 38, 239, 56, 240, 74, 241, 92, 242, 110, 243, 128, 244, 146, 245, 164, 246, 182, 247, 200};
surf_ag=news;Plane Surface(news) = {llring2,llslots};

radius =  D1s/2;
lc_circ = lcs*2;
Call Circ_ ;
llstator=newll; Line Loop(newll) = {lcirc()};


Curve Loop(newll) = {15, 1, 4, 7, 5, 8, 11, 9, -249, 219, 203, 204, 205, 210, 211, 212, 217, -248, 201, 185, 186, 187, 192, 193, 194, 199, -247, 183, 167, 168, 169, 174, 175, 176, 181, -246, 165, 149, 150, 151, 156, 157, 158, 163, -245, 147, 131, 132, 133, 138, 139, 140, 145, -244, 129, 113, 114, 115, 120, 121, 122, 127, -243, 111, 95, 96, 97, 102, 103, 104, 109, -242, 93, 77, 78, 79, 84, 85, 86, 91, -241, 75, 59, 60, 61, 66, 67, 68, 73, -240, 57, 41, 42, 43, 48, 49, 50, 55, -239, 39, 23, 24, 25, 30, 31, 32, 37, -238};
//+
surf_stator=news;Plane Surface(news) = {llstator, newll-1};


// Defining some physicals
Physical Surface("stator iron", 1000) = surf_stator;

is=0;
For k In {0:#surf_slot()-1:3}
  is=is+1;
  Physical Surface(Sprintf("slot %g",is), 100+is) = surf_slot({k:k+1});
  air_slot()+=surf_slot({k+2});
EndFor

Physical Surface("rotor shaft", 2000) = surf_shaft;
Physical Surface("rotor ring", 3000) = surf_ring;
Physical Surface("rotor ring 2", 3100) = surf_ring2;
Physical Surface("air", 3200) = {surf_ag,air_slot()};

// Some colors for aesthetics
Color SteelBlue {Surface{surf_stator, surf_shaft};}
Color SkyBlue {Surface{surf_ag,air_slot()};}
Color Magenta {Surface{surf_ring, surf_ring2};}

Color Red {Physical Surface{101:112:3};}
Color Yellow {Physical Surface{102:112:3};}
Color Green {Physical Surface{103:112:3};}
