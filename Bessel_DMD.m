clear all;

sqrt3_2=0.5*sqrt(3.0);
inv_10201=1.0/10201;
unit_image=sqrt(-1.0);

nxmax=769;
nymax=769;
totalsize=nxmax*nymax;

CenterX=(nxmax+1)/2;
CenterY=(nymax+1)/2;
Radius1=39.00;
Radius2=41.00;
Radius1_2=Radius1*Radius1;
Radius2_2=Radius2*Radius2;
Ring_Radius=ceil(Radius2)+1;

Radius3=49.00;
Radius4=51.00;
Radius3_2=Radius3*Radius3;
Radius4_2=Radius4*Radius4;
Ring_Radius2=ceil(Radius4)+1;

if Ring_Radius<Ring_Radius2
   Ring_Biggest=Ring_Radius2;
else
   Ring_Biggest=Ring_Radius;
end

Lens_Radius=0.02540/2.0;
dx=(13.68e-6);
dy=(13.68e-6);
d0=13.68e-6;
d=13.68;
size=12.68e-6;
Lambda=0.635e-6;
k0=2*pi/Lambda

initial_freq=4.65517;
final_freq=4.65517;
freq_interval=0.01;

Lattice_Lambda=13.68e-6;
M_diffraction_order=11;
N_diffraction_order=12;

tilted_angle=12*pi/180;
incident_angle=pi/7.50;
focal_length=2.50e-1;
focal_length_inv=1.0/focal_length;
focal_length2=focal_length*focal_length

p=size*(cos(tilted_angle)+1)/2;
q=size*(cos(tilted_angle)-1)/2;
gamma=sin(incident_angle)/cos(tilted_angle)-tan(tilted_angle);

Id=zeros(nxmax,nymax);

for i=1:nxmax;
    for j=1:nymax;
        dd=(i-CenterX)*(i-CenterX) + (j-CenterY)*(j-CenterY);
        if  dd<= power(Radius2,2.0) && dd>=power(Radius1,2.0)
            Id(i,j)=1;
       end

        if  dd<= power(Radius4,2.0) && dd>=power(Radius3,2.0)
            Id(i,j)=1;
       end
    end
end

for i=floor(CenterX-Ring_Biggest):ceil(CenterX+Ring_Biggest);
    i
    for j=floor(CenterY-Ring_Biggest):ceil(CenterY+Ring_Biggest);
        if i>1 && j>1 && i<nxmax && j<nymax
           if Id(i,j)==0
              if Id(i,j-1)==1 || Id(i,j+1)==1 || Id(i-1,j-1)==1 || Id(i-1,j)==1 || Id(i-1,j+1)==1 || Id(i+1,j-1)==1 || Id(i+1,j+1)==1 || Id(i+1,j)==1
                 add=0;
                 for ii=i-0.5:0.01:i+0.5
                     for jj=j-0.5:0.01:j+0.5
                         distance2=power((ii-CenterX),2.0) + power((jj-CenterY),2.0);

                         if distance2 <= Radius2_2 && distance2 >= Radius1_2
                            add=add+1;
                        end
                     end
                 end

                 if add>5100
                    Id(i,j)=1;
                end

                 add=0;
                 for ii=i-0.5:0.01:i+0.5
                     for jj=j-0.5:0.01:j+0.5
                         distance2=power((ii-CenterX),2.0) + power((jj-CenterY),2.0);

                         if distance2 <= Radius4_2 && distance2 >= Radius3_2
                            add=add+1;
                        end
                     end
                 end

                 if add>5100
                    Id(i,j)=1;
                end
             end
         else
              if Id(i,j-1)==0 || Id(i,j+1)==0 || Id(i-1,j-1)==0 || Id(i-1,j)==0 || Id(i-1,j+1)==0 || Id(i+1,j-1)==0 || Id(i+1,j+1)==0 || Id(i+1,j)==0
                 add=0;
                 for ii=i-0.5:0.01:i+0.5
                     for jj=j-0.5:0.01:j+0.5
                         distance2=power((ii-CenterX),2.0) + power((jj-CenterY),2.0);

                         if distance2 <= Radius2_2 && distance2 >= Radius1_2
                            add=add+1;
                        end
                     end
                 end

                 if add>5100
                    Id(i,j)=1;
                end

                 add=0;
                 for ii=i-0.5:0.01:i+0.5
                     for jj=j-0.5:0.01:j+0.5
                         distance2=power((ii-CenterX),2.0) + power((jj-CenterY),2.0);

                         if distance2 <= Radius4_2 && distance2 >= Radius3_2
                            add=add+1;
                        end
                     end
                 end

                 if add>5100
                    Id(i,j)=1;
                end
             end
          end
       end
    end
end

save Id Id;

Intensity=zeros(2001,2001);
cos_thita=cos(2*tilted_angle);
k_costhitain_f=k0*d0*sin(2*tilted_angle)/focal_length*(1.0e-6);
k_sinthitaIn=k0*sin(2*tilted_angle)
k0_focal_length=k0/focal_length
k0_time_r0=k0*(Radius1+Radius2)/2*d0
k0_time_r1=k0*(Radius3+Radius4)/2*d0
k0_Lens_Radius_focal_length=k0*Lens_Radius/focal_length*(1.0e-06)
Factor1=((Radius2-Radius1)*d0)
Factor2=unit_image*sin(k0_focal_length*(Radius1+Radius2)/2*d0*0.5*Factor1)
Factor3=cos(k0_focal_length*(Radius1+Radius2)/2*d0*0.5*Factor1)
Factor4=1.0/sin(2*tilted_angle)
Factor5=1.0/focal_length/focal_length
Factor6=sin(2*tilted_angle)/focal_length
Factor7=k0*0.5*Factor1
sinthitain2=sin(2*tilted_angle)*sin(2*tilted_angle)
delta_r0=((Radius2-Radius1)*d0)
delta_r1=((Radius4-Radius3)*d0)
r0=(Radius1+Radius2)/2*d0
r1=(Radius3+Radius4)/2*d0
r0_minus_delta_r=(Radius1+Radius2)/2*d0-0.5*Factor1
r0_plus_delta_r=(Radius1+Radius2)/2*d0+0.5*Factor1
delta_Lens_Radius=Lens_Radius/6000
delta_phi=pi/900.0
delta_Lens_Radius2=delta_Lens_Radius*delta_Lens_Radius
Sinc_alpha=pi*size/Lambda
k0_size=k0*(0.5*size)
beta=0.0;

Amplitude_Lens=zeros(5001,5001);
Id2=zeros(5001,5001);
for i=203:1365;
    i2=(i-1001)*(i-1001);
    i_size=(i-1001)*(delta_Lens_Radius)

    for j=420:1582;
        r=sqrt(i2+(j-1001)*(j-1001));
        j_size=(j-1001)*(delta_Lens_Radius);

        if r>0
           beta=acos((i-1001)/r);
        end

        if r*delta_Lens_Radius<=Lens_Radius
           Id2(i,j)=1;
       end

        %%accumulation=0;
        Amplitude=0.0;
        for i0=floor(CenterX-Ring_Biggest):ceil(CenterX+Ring_Biggest)
            i02=(i0-CenterX)*(i0-CenterX);
            i12=((i0-CenterX)*d0-i_size);
            i12_2=i12*i12;


            for j0=floor(CenterY-Ring_Biggest):ceil(CenterY+Ring_Biggest)
                if Id(i0,j0)==1
                   Radius2=(i02+(j0-CenterY)*(j0-CenterY));
                   Radius=sqrt(Radius2);
                   cos_thita=(i0-CenterX)/Radius;

                   j12=((j0-CenterY)*d0-j_size);
                   j12_2=j12*j12;

                   s=sqrt(focal_length2+i12_2+j12_2);
                   s_inv=1.0/s;

                   sin_thita_x=i12*s_inv;
                   sin_thita_y=j12*s_inv;

                   cos_R_f=focal_length*s_inv;
                   thita_R_f=acos(cos_R_f);

                   %%sin_ns=sqrt(i12_2+j12*j12)*focal_length_inv;
                   %%Radius_Dot_r=((i0-CenterX)*(i-1)+(j0-CenterY)*(j-1));
                   %%qy=abs(j12)*k0_focal_length;

                   Amplitude=Amplitude+exp(unit_image*(k_sinthitaIn*Radius*d0*cos_thita))*...
                                       exp(unit_image*k0*s)*...
                                       sinc(k0_size*sin_thita_x)*sinc(k0_size*sin_thita_y)*...
                                       s_inv;
               end
            end
        end

        Amplitude_Lens(i,j)=Amplitude;
    end
end

save Amplitude_Lens Amplitude_Lens;
%%load Amplitude_Lens;
for i=1:2001;  %%500
    i2=(i-1001)*(i-1001);
    i_size=(i-1001)*(1.0e-6)

    for j=1:2001;
        r=sqrt(i2+(j-1001)*(j-1001));
        j_size=(j-1001)*(1.0e-6)

        if r>0
           beta=acos((i-1001)/r);
       end

        Amplitude=0.0;
        for i0=203:1365
            i02=(i0-1001)*(i0-1001);

            for j0=420:1582;
                if Id2(i0,j0)==1
                   j02=(j0-1001)*(j0-1001);
                   Lens_r2=(i02+j02)*delta_Lens_Radius2;
                   Lens_r=sqrt(Lens_r2);

                   Amplitude=Amplitude+focal_length_inv*...
                                       (Amplitude_Lens(i0,j0))*...
                                        exp(unit_image*0.5*k0_focal_length*Lens_r2)*...
                                        besselj(0,k0_focal_length*Lens_r*r*(1.0e-6))*...
                                        Lens_r*delta_Lens_Radius;
               end
            end
        end

        Intensity(i,j)=(Amplitude);
    end
end

save Intensity Intensity;

figure(1)
color_low=0.0;
color_top=max(max(Amplitude_Lens(:,:).*conj(Amplitude_Lens(:,:))));
clf;
pcolor(abs(Amplitude_Lens.*Amplitude_Lens));
axis('equal');
axis([1 2001 1 2001]);
caxis([color_low color_top]);
shading interp;
colorbar;
xlabel('x-axis');
ylabel('y-axis');

figure(2)
color_low=0.0;
color_top=max(max(Intensity(:,:).*conj(Intensity(:,:))));
clf;
pcolor(abs(Intensity.*Intensity));
axis('equal');
axis([1 2001 1 2001]);
caxis([color_low color_top]);
shading interp;
colorbar;
xlabel('x-axis');
ylabel('y-axis');

save Intensity Intensity;

figure(3)
i=1:2001;
plot(i,Amplitude_Lens(784,i).*conj(Amplitude_Lens(784,i))./(max(Amplitude_Lens(784,:).*conj(Amplitude_Lens(784,:)))),'b');
hold on;
xlabel('Position (um)');
ylabel('Normalized Intensity');

figure(4)
i=1:5000;
plot(i,Intensity(i,1001).*conj(Intensity(i,1001))./(max(IntensityAll(:,1001).*conj(IntensityAll(:,1001)))),'b');
hold on;
xlabel('Position (um)');
ylabel('Normalized Intensity');

