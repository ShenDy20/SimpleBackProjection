tic;
SPACE=zeros(512,512,600); % the space where the source is located
for I=1:length(ZG)
    v1=[XG(I,2)-XG(I,1) YG(I,2)-YG(I,1) -ZG(I)/55];
    co=1-511*(1/EG(I,2)-1/(EG(I,1)+EG(I,2)));
    l=norm(v1);
    q=XG(I,2)-XG(I,1);
    w=YG(I,2)-YG(I,1);
    e=-ZG(I)/55;
    for a=1:512
        for b=1:512 % for each (x,y), a certain cone only passes through certain Z 
            a1=(XG(I,1)+128-a);  
            b1=(YG(I,1)+128-b);  % add the x y offset of the detector
            A=e^2-(co^2)*(l^2);
            B=e*(a1*q+b1*w);
            C=(a1*q+b1*w)^2-(co^2)*(l^2)*(a1^2+b1^2);
            z1=-(-B+sqrt(B^2-4*A*C))/(2*A);
            z2=-(-B-sqrt(B^2-4*A*C))/(2*A);
            if(z1>=1&&z1<=600&&isreal(z1))
                v2=[a1 b1 -z1];
                if dot(v1,v2)*co>0
                SPACE(a,b,round(z1))=1+ SPACE(a,b,round(z1));
                end
            end
            if(z2>=1&&z2<=600&&isreal(z2)&&(round(z2)~=round(z1))) 
                v3=[a1 b1 -z2];
               if dot(v1,v3)*co>0
                SPACE(a,b,round(z2))=1+ SPACE(a,b,round(z2));
               end
            end
        end
    end
    
end
fx=[];
fy=[];
fz=[];
gx=[];
gy=[];
gz=[];
hx=[];
hy=[];
hz=[];
for c=1:600
    for a=1:512
        for b=1:512
             if SPACE(a,b,c)>10  %  threshold for determining whether there's a source in certain pixel
                hx=[hx;a];
                hy=[hy;b];
                hz=[hz;c];
            end
        end
    end
end

toc;
 
  save ('u.txt','hx','-ascii')
 save ('v.txt','hy','-ascii')
 save ('w.txt','hz','-ascii')

                