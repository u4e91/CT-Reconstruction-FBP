%丁昊妍 520021910379

function [Img] = fbp(projection,filter_choose, sensor_length,pixelsize, app)
    detector_size = sensor_length;
    M = size(projection,1);
    N_d = M;
    TT = linspace(-N_d*detector_size/2,N_d*detector_size/2,N_d+1);    %探测器边界坐标
    for i=1:N_d
        DD(i) = (TT(i)+TT(i+1))/2;    %#ok<NASGU,AGROW> %探测器通道的中点坐标
    end
    
    Img = zeros(M);
    [ImgRow,ImgCol] = size(Img);
    Detectorsize = detector_size;
    
    if filter_choose == 1
        %和滤波器卷积
        dS=1;
        rayNum = sensor_length;
        ray_x=-(rayNum-1)/2 * dS:dS:(rayNum-1)/2 * dS;
        ray_x=ray_x';

        l=length(ray_x);

        N=M;
        t = linspace(-N/2,N/2-1,N);%生成线性间距向量
        filt = 0.0085*(sinc(t)/2-sinc(t/2).^2/4);%得到滤波的卷积核
        
        p = projection;
        projection0 = zeros(l+2*N,180);
        projection0(N+1:l+N,:) = p;
        for n = 1:180
            temp(:,n) = conv(projection0(:,n),filt); %#ok<AGROW>
        end
        projection0=temp(round(3*M/2:3*M/2+l),:);%投影滤波
        projection0_show = imresize(projection0, [256,256]);
        imshow(projection0_show,[],'Parent',app.UIAxes4)
        xlabel('\theta(degrees)','Parent',app.UIAxes4);
        ylabel('x','Parent',app.UIAxes4);
        drawnow
        %卷积结束

        %Img = zeros(M);
        %Detectorsize = detector_size;
        %[ImgRow,ImgCol] = size(Img);
        width = 2^nextpow2(size(projection0,1))*16;
        pro_fft = fft(projection0,width);

        %R_L = (1/Detectorsize)*[0:(width/2), width/2-1:-1:1]'/width; %1/2tau,tau为采样间距
        filter_Hamming = 0.54-0.46*cos([-width/2+1:0,1:width/2]'.*(2*pi/(width)));
        pro_filter = zeros(width,size(projection0,2));
        for i = 1:1:size(projection0,2)
            pro_filter(:,i)=pro_fft(:,i).*filter_Hamming;
            %pro_filter(:,i) = pro_fft(:,i).*R_L.*filter_Hamming;
        end
        pro_ifft = real(ifft(pro_filter));
    else
        pro_ifft = projection;
    end
    
    
    for i = 1:1:size(projection,2)
        for c = 1:ImgCol
            for r = 1:ImgRow
                x = (c-(ImgCol+1)/2).*pixelsize(1);
                y = ((ImgRow+1)/2-r).*pixelsize(2);
                s = x*cos((i-1)*pi/180)+y*sin((i-1)*pi/180);
                s = s/Detectorsize + N_d/2 + 1;
                n = floor(s);
                if n > 0 && n <= N_d
                    q = pro_ifft(n,i);
                    Img(r,c) = Img(r,c)+q;
                end
            end
        end
        str = ['重建进度：',num2str(round(i/size(projection,2)*100)),'%'];
        app.updateGUI(str);
        imshow(Img,[],'Parent',app.UIAxes3);
        drawnow
    end
    Img = Img*pi/size(projection,2);
    
end

