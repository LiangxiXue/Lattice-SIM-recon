
%Reconstruction of super-resolution imaging:
tic
close all; clear all; clc; 
Path_name='G:\Experiment\Lattice_SIM_test_Data\'; %Path name; 

Pixel=19.5; %Pixel size (nm); 
Wave=532e-9; %Wavelength of the light(m); 
k=2*pi/Wave; %Wave vector of red light;  
Apodization=0; %"1" to perform apodization, '0' not to perform apodization; 
m_s=0.7; %modulation factor; 
m_t=0.7; %modulation factor;
Image_wf=zeros(1000,1000); 
Image_SIM=zeros(1000,1000); 

for Num_y=2:1:2,
   for Num_x=2.5:2.5,
           y_start=1000*(Num_y-1)+2; y_end=1000*(Num_y-0)+25; if Num_y==3,y_end=1000*(Num_y-0); end
           x_start=1000*(Num_x-1)+3; x_end=1000*(Num_x-0)+26;  if Num_x==4,x_end=1000*(Num_x-0); end
           
           Image_raw=zeros(y_end-y_start+1,x_end-x_start+1,1); %Innitialization;

            %Load image 1:
            for ii=1:1:5,
                File_name=strcat(num2str(ii),'.tiff'); %Filename;
                TifLink_SIM_ud= Tiff([Path_name,File_name], 'r'); 
                Pic=TifLink_SIM_ud.read();     Pic=double(Pic); %read raw images;
                Image_raw(:,:,ii)=Pic(y_start:y_end,x_start:x_end); %Raw image;
                %mean(mean(Image_raw))
                %figure(1); imagesc(squeeze(Image_raw(:,:,ii))); colormap('gray');
                %pause(0.5);
                
            end

            %Intensity normalization:
            Intensity_ave=mean(mean(mean(Image_raw))); 
            for ii=1:1:5,
                Pic=squeeze(Image_raw(:,:,ii)); 
                Pic=Pic/mean(mean(Pic))*Intensity_ave; 
                Image_raw(:,:,ii)=Pic; %Normalized intensity; 
                %figure(1); imagesc(Pic); colormap('gray');      
            end

            %3x2 phase-shifting-reconstruction:
            Coef_matrix=zeros(2,2);
            Num=1; 
            for ii=1:1:2,  %(y);
                Pha_y=2*pi/3*(ii-1); %Phase-shifts in the y direction; 
                for jj=1:1:3, %(x);     
                    Pha_x=2*pi/3*(jj-1); %Phase-shifts in the x direction; 
                    Coef_matrix(Num,1)=2;  
                    Coef_matrix(Num,2)=exp(1i*Pha_x); Coef_matrix(Num,3)=exp(-1i*Pha_x); %Coefficient; 
                    Coef_matrix(Num,4)=exp(1i*Pha_y); Coef_matrix(Num,5)=exp(-1i*Pha_y); %Coefficient; 
                    Num=Num+1; 
                end
            end
            Coef_matrix(5,4)=exp(1i*2*pi/3*2); Coef_matrix(5,5)=exp(-1i*2*pi/3*2); 
            Coef_LS=Coef_matrix(1:5,:); %Coefficient matrix; 

            %figure(5); imagesc(squeeze(Image_raw(:,:,3))); colormap('gray');

            Image_fre=zeros(size(Image_raw));  % Images along different diffraction orders of the structured illumination; 
            Coef_Inv=inv(Coef_LS); %inverse matrix; 

            for ii=1:1:5, 
                for jj=1:1:5, 
                     Image_fre(:,:,ii)=Image_fre(:,:,ii)+Coef_Inv(ii,jj)*Image_raw(:,:,jj); 
                end
            end
            %figure(5); imagesc(sum(Image_raw(:,:,1:1:5),3)); colormap('gray');

            %The components along different diffraction orders:
            fAo= squeeze(Image_fre(:,:,1))/2; %The dc term;
            fAp= squeeze(Image_fre(:,:,2));  %The +1st order in the s direction; 
            fAm=squeeze(Image_fre(:,:,3));  %The -1st order in the s direction; 
            fBo= squeeze(Image_fre(:,:,1))/2; %The dc term; 
            fBp= squeeze(Image_fre(:,:,4)); %The +1st order in the t direction; 
            fBm=squeeze(Image_fre(:,:,5)); %The -1st order in the t direction; 
            
            figure(1); imagesc(abs(fAo));          
            figure(2); imagesc(abs( fAp));
            factor1=mean(mean(abs(fAm)./abs(fAo)));
            factor2=mean(mean(abs(fBp)./abs(fBo)));

            %Coordinates in spatial space: 
            [y_length,x_length]=size(fAo); %The size of hologram; 
            x=linspace(-x_length/2,x_length/2-1,x_length);  %x coordinate(micro); 
            y=linspace(-y_length/2,y_length/2-1,y_length); %y coordinate(micro); 
            [xx,yy]=meshgrid(x,y); %Two-dimensional coordinates; 
            
            %Coordinates in Fourier space: 
            u=linspace(-x_length/2,x_length/2-1,x_length)/x_length;  u=u*2*pi; %Spectrum coordinate; 
            v=linspace(-y_length/2,y_length/2-1,y_length)/y_length; v=v*2*pi; %Range: 0~2pi; 
            [uu,vv]=meshgrid(u,v); %Generate the two-dimensional coordinates; 
            
            %%spectra are moved to their correct position in Fourier space, added up, and compensated for the frequency damping 
            %Deconvolution through Wiener-filtering:
            Deno1=0; Num1=0;
            
                          %Spectrum of the (st,0) term:
                            [w0_x,w0_y]=Single_Carrier_frequency_detection(fAp);
                            Mask0=double(sqrt((uu+0).^2+(vv+0).^2)<sqrt(w0_x^2+w0_y^2)/1)*1/0.85; 
                            OTF=exp(-((uu+0).^2+(vv+0).^2)/(w0_x^2+w0_y^2)); %the OTF of the imaging system; 
                            Freq_0= fftshift(fft2(fftshift(fAo+fBo)));%.*Mask0;  
                            Deno1=Deno1+OTF.^2/1; %
                            Num1=Num1+Freq_0.*conj(OTF); %

                          %Spectrum of the (s,1) term;
                            [w0_x,w0_y]=Single_Carrier_frequency_detection(fAp);
                            Mask1=double(sqrt((uu+w0_x/1).^2+(vv+w0_y/1).^2)<sqrt(w0_x^2+w0_y^2)/0.85); 
                            OTF=exp(-((uu+w0_x).^2+(vv+w0_y).^2)/(w0_x^2+w0_y^2)); %the OTF of the imaging system;                       
                            Sim_R=exp(-1i*(w0_x*xx+w0_y*yy)); %Reference wave; 
                            fAp=fAp.*Sim_R;           
                            Freq_Ap= fftshift(fft2(fftshift(fAp)));             
                            Freq_Ap=Freq_Ap.*exp(-1i*angle(Freq_Ap(513,513)./Freq_0(513,513)));%.*Mask1;            
                            Deno1=Deno1+OTF.^2; %
                            Num1=Num1+Freq_Ap.*conj(OTF)/m_s;

                            %Spectrum of the (s,-1) term;
                            [w0_x,w0_y]=Single_Carrier_frequency_detection(fAm); 
                            Mask2=double(sqrt((uu+w0_x/1).^2+(vv+w0_y/1).^2)<sqrt(w0_x^2+w0_y^2)/0.85); 
                            OTF=exp(-((uu+w0_x).^2+(vv+w0_y).^2)/(w0_x^2+w0_y^2)); %the OTF of the imaging system; 
                            Sim_R=exp(-1i*(w0_x*xx+w0_y*yy)); %Reference wave; 
                            fAm=fAm.*Sim_R; 
                            Freq_Am= fftshift(fft2(fftshift(fAm))); 
                            Freq_Am=Freq_Am*exp(-1i*angle(Freq_Am(513,513)./Freq_0(513,513)));%.*Mask2;
                            Deno1=Deno1+OTF.^2; %
                            Num1=Num1+Freq_Am.*conj(OTF)/m_s;

                            %Spectrum of the (t,1) term;
                            [w0_x,w0_y]=Single_Carrier_frequency_detection(fBp);
                            Mask3=double(sqrt((uu+w0_x/1).^2+(vv+w0_y/1).^2)<sqrt(w0_x^2+w0_y^2)/0.85); 
                            OTF=exp(-((uu+w0_x).^2+(vv+w0_y).^2)/(w0_x^2+w0_y^2)); %the OTF of the imaging system; 
                            Sim_R=exp(-1i*(w0_x*xx+w0_y*yy)); %Reference wave; 
                            fBp=fBp.*Sim_R; 
                            Freq_Bp= fftshift(fft2(fftshift(fBp))); 
                            Freq_Bp=Freq_Bp*exp(-1i*angle(Freq_Bp(513,513)./Freq_0(513,513)));%.*Mask3;
                            Deno1=Deno1+OTF.^2; %
                            Num1=Num1+Freq_Bp.*conj(OTF)/m_t;

                            %Spectrum of the (t,-1) term;
                            [w0_x,w0_y]=Single_Carrier_frequency_detection(fBm); %%%GRATE                          
                            Mask4=double(sqrt((uu+w0_x/1).^2+(vv+w0_y/1).^2)<sqrt(w0_x^2+w0_y^2)/0.85); 
                            OTF=exp(-((uu+w0_x).^2+(vv+w0_y).^2)/(w0_x^2+w0_y^2)); %the OTF of the imaging system; 
                            Sim_R=exp(-1i*(w0_x*xx+w0_y*yy)); %Reference wave; 
                            fBm=fBm.*Sim_R; 
                            Freq_Bm= fftshift(fft2(fftshift(fBm))); 
                            Freq_Bm=Freq_Bm*exp(-1i*angle(Freq_Bm(513,513)./Freq_0(513,513)));%.*Mask4;
                            Deno1=Deno1+OTF.^2; %
                            Num1=Num1+Freq_Bm.*conj(OTF)/m_t;
            
            %figure(100); imagesc(log(1+abs(Num1)));
            %figure(200); imagesc(log(1+abs(Deno1))); 
           
            Freq_total=Num1./(Deno1+0.04); %0.04, the omega dampens the degree of compensation especially in regions with small values in the denominator
            
            %Apodization (Optics Express 21 2032-2049 2013) of the combined spectrum to compensate for ringing artifact:
            ApoFunc=1;
            if Apodization==1
                    ApoMask=double(Mask0+Mask1+Mask2+Mask3+Mask4<1);  %Reversed spectrum domain;
                    DistApoMask = bwdist(ApoMask); 
                    maxApoMask = max(max(DistApoMask)); 
                    ApoFunc = double(DistApoMask./maxApoMask).^0.4; 
            end
            %figure(100); imagesc(ApoFunc);
            
            Freq_total=Freq_total.*ApoFunc; %Apodized SIM spectrum; 
            
            %Zero-filling in the spectrum in order to reconstruct an image with double size
            Freq_total_DS=zeros(2*y_length,2*x_length); %Double size;
            Freq_total_DS(0.5*y_length+1:1.5*y_length,0.5*x_length+1:1.5*x_length)=Freq_total;
            Freq_0_DS=zeros(2*y_length,2*x_length); %Double size;
            Freq_0_DS(0.5*y_length+1:1.5*y_length,0.5*x_length+1:1.5*x_length)=Freq_0;
            
            %Reconstructed wide-field and SIM image
            I_wf=abs(fftshift(ifft2(fftshift(Freq_0_DS)))); %Reconstructed SIM image; 
            I_SIM=abs(fftshift(ifft2(fftshift(Freq_total_DS)))); %Reconstructed SIM image; 
            
            %figure(100); imagesc(abs(Freq_Ap.*Mask1)); caxis([0,20000]); colormap('hot');
            %figure(200); imagesc(abs(Freq_total.*ApoFunc)); caxis([0,40000]); colormap('hot')
            figure(400); imagesc(I_wf); colormap('gray');  caxis([0,50]); colorbar; title('wide-field image'); colormap('hot');
            figure(500); imagesc(I_SIM); colormap('gray');  caxis([0,70]); colorbar; title('Reconstructed SIM image'); colormap('hot');
                       

%             Image_wf(y_start:y_end,x_start:x_end)=real(fAo); % Raw image;
%             Image_SIM(y_start:y_end,x_start:x_end)=abs(Dsum); % SIM image;
            
   end
end

toc
%Write the reconstruction in TIFF image: 
I_SIM=I_SIM/40;
cd (Path_name); imwrite(I_wf/20*1,'Wide-field.tiff','tiff'); 
cd (Path_name); imwrite(I_SIM,'SIM.tiff','tiff'); 

% I_ref=I_SIM; 
% cd G:\Experiment\2022_06_08_Lattice\Lattice\2\2; load I_ref I_ref; 
% figure(600); imagesc(abs(I_ref(500:600,500:600))); 
% figure(700); imagesc(abs(I_SIM(500:600,500:600))); 
% 
% SSIM_val = ssim(I_SIM, I_ref), 
  