%Purpose:get the carrier frequency from the carrier interferogram;
%input:  'Holo',the intensity distribution of carrier-frequency interferogram;
%Output: 'w0x1','w0y1', Carrier frequency in x or y direction; Unit:rad/pixel;
         %'w0x2','w0y2', for the second carrier frequency;
       
function [w0x,w0y]=Single_Carrier_frequency_detection(Holo)

       %Holo=imresize(Holo,[3072,3840],'cubic'); %resize; 
       %Holo=Holo(1:1016,1:1272); %resize; 
       Carr_ratio=1/15; 
       
      %Zero padding: 
        [y_length,x_length]=size(Holo); %Get the size of hologram; 
        Fre_I=fftshift(fft2( fftshift(Holo) )); %Frequency spectrum; 
        Fre_I(Fre_I~=Fre_I)=1; %Eleminate the null data; 
        Fre_I=abs(Fre_I); %Filtering; 
        
    %Coordinates: 
        u=linspace(-x_length/2,x_length/2-1,x_length)/x_length;  u=u*2*pi; %Spectrum coordinate; 
        v=linspace(-y_length/2,y_length/2-1,y_length)/y_length; v=v*2*pi; %Range: 0~2pi; 
        [uu,vv]=meshgrid(u,v); %Generate the two-dimensional coordinates; 
        Angle_pol=angle(uu+1i*vv); %Polar angle; from -pi ~ pi; 
        
        Mask=double( sqrt(uu.^2+vv.^2)>=max(u)*Carr_ratio );  Fre_I=Fre_I.*Mask; %Filtering the zeroth order; 
        %figure(1000); imagesc(u,v,Angle_pol); grid on; 
        
    %For one carrier-frequenccy: (Frequency spectrum in quadrature 1): 
        Mask1=(Angle_pol>-2*pi & Angle_pol<=2*pi); 
        uu1=uu(Mask1); vv1=vv(Mask1); % Quadrature 1; 
        Fre_sel_x=Fre_I(Mask1);  Fre_sel_y=Fre_I(Mask1); % Quadrature 1; 
        w0x=uu1(Fre_sel_x==max(max(Fre_sel_x))); %Carrier frequency in x direction; Unit:rad/pixel; 
        w0y=vv1(Fre_sel_y==max(max(Fre_sel_y))); %Carrier frequency in y direction; Unit:rad/pixel; 
        
       %figure(10); imagesc(u,v,log(Fre_I)); grid on; %Display the frequency spectrum;   
        
end