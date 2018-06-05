function cond = ORdConductances(celltype)
%endo = 0, epi = 1, M = 2

GNa = 75;

GNaL=0.0075;
if celltype==1
    GNaL=GNaL*0.6;
end

Gto=0.02;
if celltype==1
    Gto=Gto*4.0;
elseif celltype==2
    Gto=Gto*4.0;
end

PCa=0.0001;
if celltype==1
    PCa=PCa*1.2;
elseif celltype==2
    PCa=PCa*2.5;
end

GKr=0.046;
if celltype==1
    GKr=GKr*1.3;
elseif celltype==2
    GKr=GKr*0.8;
end

GKs=0.0034;
if celltype==1
    GKs=GKs*1.4;
end

GK1=0.1908;
if celltype==1
    GK1=GK1*1.2;
elseif celltype==2
    GK1=GK1*1.3;
end

Gncx=0.0008;
if celltype==1
    Gncx=Gncx*1.1;
elseif celltype==2
    Gncx=Gncx*1.4;
end

Pnak=30;
if celltype==1
    Pnak=Pnak*0.9;
elseif celltype==2
    Pnak=Pnak*0.7;
end

GKb=0.003;
if celltype==1
    GKb=GKb*0.6;
end

PNab=3.75e-10;

PCab=2.5e-8;

GpCa=0.0005;

Pjrel = 1;

Pjup = 1;

cond = [GNa , GNaL , Gto , PCa , GKr , GKs , GK1 , Gncx , Pnak, GKb, PNab, PCab, GpCa, Pjrel, Pjup];