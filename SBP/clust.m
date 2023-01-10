tic;
XG=[];
 YG=[];
 TG=[];
 EG=[];
 ZG=[];
parfor HH=0:30
    if HH<10
        
filename=['./atoms00' num2str(HH) '.txt'] ;
    else
        filename=['./atoms0' num2str(HH) '.txt'] ;
    end
[x,y,t,e]=textread(filename,'%n %n %f %f','delimiter', ' ');


for GEN=0:26
X=zeros(10000,1);
Y=zeros(10000,1);
T=zeros(10000,1);
E=zeros(10000,1);
    for c=1:10000
        cons=GEN*10000+c;
        X(c)=x(cons);
        Y(c)=y(cons);
        T(c)=t(cons); 
        E(c)=e(cons);
    end
CL=zeros(40,5000);
TL=CL;
EL=CL;
XL=CL;
YL=CL;
start=0;
labels=dbscan(T,100,2);     % time clustering
numGroups = length(unique(labels)); 
pos=zeros(numGroups-1); 
    for I=1:10000
        if labels(I)~=-1
       pos(labels(I))=pos(labels(I))+1;
       CL(pos(labels(I)),labels(I))=I;
       TL(pos(labels(I)),labels(I))=T(I);
       EL(pos(labels(I)),labels(I))=E(I);
       XL(pos(labels(I)),labels(I))=X(I);
       YL(pos(labels(I)),labels(I))=Y(I);
        end
    end
   
  for J=1:numGroups-1
     k=1;
     tempx=[];
     tempyy=[];
     tempt=[];
     tempe=[];
     temp=[];
     while XL(k,J)~=0
         tempx=[tempx XL(k,J)];
         tempyy=[tempyy YL(k,J)];
         tempt=[tempt TL(k,J)];
         tempe=[tempe EL(k,J)];
         temp=[temp; XL(k,J)+YL(k,J)];
         k=k+1;
         if k==40
             break; 
         end
     end
     temp=sort(temp);
     ps=dbscan(temp,1.1,1);  % space clustering
     group=length(unique(ps));
     if group==1
         continue
     end
     p=zeros(group+2);
     cl=zeros(10,10);
     tl=cl;
     el=cl;
     xl=cl;
     yl=cl;
     for m=1:k-1
         tem=ps(m)+2;
         pos(tem)=pos(tem)+1;
         n=pos(tem);
         cl(n,tem)=m;
         tl(n,tem)=tempt(m);
         el(n,tem)=tempe(m);
         xl(n,tem)=tempx(m);
         yl(n,tem)=tempyy(m);
     end
     
      if group==2  % reserve compton incidents with only two energy deposit
       xc1=xl(:,3);
       xc2=xl(:,4);
       yc1=yl(:,3);
       yc2=yl(:,4);
       yy=[mean(yc1(yc1~=0)) mean(yc2(yc2~=0))];
       xx=[mean(xc1(xc1~=0)) mean(xc2(xc2~=0))];
       tc1=tl(:,3);
       tc2=tl(:,4);
       tt=[min(tc1(tc1~=0)) min(tc2(tc2~=0))];
       ee=[sum(el(:,3)) sum(el(:,4))];
       ee=sort(ee);
       tt=sort(tt);

 if ee(1)+ee(2)>50&&ee(1)+ee(2)<70&&abs(1-511*(1/ee(2)-1/(ee(1)+ee(2))))<1 
     % make sure the sum of two energy deposit is close to the feature
     % energy of the source, so that the scattering only happens once
        if sum(el(:,3))*2>ee(1)+ee(2)
        yyy=yy(1);
        yy(1)=yy(2);
        yy(2)=yyy;
        xxx=xx(1);
        xx(1)=xx(2);
        xx(2)=xxx;
       
        end
        % XG YG ZG EG save the location and energy information of compton incidents 
        YG=[YG;yy];
        TG=[TG;tt];
        EG=[EG;ee];
        XG=[XG;xx];
        ZG=[ZG;(tt(2)-tt(1))*18.4837]; % the number 18.4837 is derived from experiment as the drifting speed
         start=start+1;
       end
     end
      
  end


end
end
toc;