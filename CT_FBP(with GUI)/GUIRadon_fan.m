%丁昊妍 520021910379

imshow([],'Parent',app.UIAxes2);
imshow([],'Parent',app.UIAxes3);
imshow([],'Parent',app.UIAxes4);
xlabel([],'Parent',app.UIAxes2);
ylabel([],'Parent',app.UIAxes2);
xlabel([],'Parent',app.UIAxes4);
ylabel([],'Parent',app.UIAxes4);

img = app.LoadButton.UserData;
filter_choose = app.filter_choose;
if app.P_choose == 1
    M=256;
    [img,E] = phantom('Modified Shepp-Logan', M); % 创建一个 shepp-logan 模型，原图像
end  

beta=0:1:359;%射线源与探测器旋转角度，默认逆时针旋转
[M,~]=size(img);
dS=1;%探测器在物体中心线处的间隔（探测器实际上间隔会大于dS）
rayNum=ceil(sqrt(2)*M);%探测器数=射线数
distance=2*M;%射线源到物体中心的距离，SOD

R=Projection(img,beta,distance,rayNum,dS, app);%正投影，因为用了interp2比较慢
R1 = imresize(R, [256,256]);
imshow(R1,[],'Parent',app.UIAxes2);
xlabel('\theta(degrees)','Parent',app.UIAxes2);
ylabel('x','Parent',app.UIAxes2);
drawnow

[recon]=FanBeamFBP(M,M,beta,R,rayNum,dS,distance,filter_choose, app);%滤波反投影
imshow(recon,[],'Parent',app.UIAxes3);



function u=U(D,x,y,beta)%公式中的U
    u=(D+x*sin(beta)-y*cos(beta))/D;
end

function proj = Projection(img, beta,distance,rayNum,dS, app)%正投影过程
    proj = zeros(rayNum, length(beta));
    [M,~]=size(img);
    M2=(M-1)/2;
    [X, Y] = meshgrid(-M2 : 1 : M2);%图片坐标系重定位，以图像中心为(0,0)
    ray_x=-(rayNum-1)/2 * dS:dS:(rayNum-1)/2 * dS;%各光线在物体中心坐标系中与射线源中心线垂直线上的各值，dS为各探测器在该线上的间隔
    CY=-M2*sqrt(2):1:M2*sqrt(2);%垂直x轴线采样点的y,图像为正方形且采样间隔为1，因此设采样点数为根号2*M+1
    CX=zeros([1,length(CY)]);%垂直x轴线采样点的x
    for i=beta%360角全投影,默认逆时针旋转
        proj_1=zeros(rayNum,1);%一个角度的投影
        r=deg2rad(i);%角度转弧度
        CX_1=CX.*cos(-r)+CY.*sin(-r);
        CY_1=CY.*cos(-r)-CX.*sin(-r);%射线源旋转i度后采样点跟着旋转
        for j=1:length(ray_x)%每个角度有ray_x根线
            gamma = atan(ray_x(j)/distance);%gamma角，即射线与穿过源射线的夹角
            theta = r+gamma;%cos(theta)*x+sin(theta)*y=t中的theta
            t = (distance*ray_x(j))/sqrt(distance*distance+ray_x(j)*ray_x(j));%cos(theta)*x+sin(theta)*y=t中的t
            bias_x = t*cos(theta);%扇形束采样点的x偏置
            bias_y = t*sin(theta);%扇形束采样点的y偏置
            CX_real = CX_1.*cos(-gamma)+CY_1.*sin(-gamma)+bias_x;
            CY_real = CY_1.*cos(-gamma)-CX_1.*sin(-gamma)+bias_y;%采样点旋转gamma角并加上偏置得到射线采样点
            inter = interp2(X, Y, img, CX_real, CY_real, 'linear',0)';%线性插值得到各采样点的值,如果可以最好不要使用interp函数，很慢
            proj_1(j) = sum(inter,1); %#ok<UDIM>
        end
        proj(:,i+1)=proj_1;        
        str = ['扫描进度：',num2str(round(i/359*100)),'%'];
        app.updateGUI(str);
        drawnow
        
        R1 = imresize(proj, [256,256]);
        imshow(R1,[],'Parent',app.UIAxes2);
        xlabel('\theta(degrees)','Parent',app.UIAxes2);
        ylabel('x','Parent',app.UIAxes2);
        
    end
end

function [recon] = FanBeamFBP(sizeM,sizeN,beta,R,rayNum,dS,distance,filter_choose, app)%反投影过程
    M=sizeM;%重建图像的行数
    N=sizeN;%重建图像的列数
    ray_x=-(rayNum-1)/2 * dS:dS:(rayNum-1)/2 * dS;
    ray_x=ray_x';
    weight_d=ones(rayNum,1).*distance;
    weight=weight_d./(sqrt(weight_d.*weight_d+ray_x.*ray_x));%得到权重
    clear weight_d;
    l=length(ray_x);
    clear ray_x;
    for i = 1:length(beta)
        R(:,i)=R(:,i).*weight;%投影加权
    end
    clear weight;
    if filter_choose == 1
        t=linspace(-N/2,N/2-1,N);
        filt=0.0085*(sinc(t)/2-sinc(t/2).^2/4);%得到滤波的卷积核
        p=R;
        R = zeros(l+2*N,360);
        R(N+1:l+N,:) = p;
        for n = 1:360
            temp(:,n) = conv(R(:,n),filt); %#ok<AGROW>
        end
        R=temp(3*M/2:3*M/2+l,:);%投影滤波
        R_filt1 = imresize(R, [256,256]);
        imshow(R_filt1,[],'Parent',app.UIAxes4);
        xlabel('\theta(degrees)','Parent',app.UIAxes4);
        ylabel('x','Parent',app.UIAxes4);
        drawnow
        
    end
    clear p;
    recon=zeros([M,N]);%重建图像
    M2=(M-1)/2;
    dbeta=deg2rad(1);
    for i=beta%射线源旋转每一个角度
        r=deg2rad(i);%角度转弧度
        for x=1:M
            for y=1:N%对于重建图像中的每一个点
                sp=distance*((y-M2-1)*cos(r)+(x-M2-1)*sin(r))/(distance+(y-M2-1)*sin(r)-(x-M2-1)*cos(r));
                %得到s'，这里的推导可以试着自己推一推，过(x,y)和射线源的点和D‘2(该线随β变化而变化)的交点与O连接线段的长度。再注意：这段循环中的x和y并不是正常的x和y，x表示行数而y表示列数，与我们习惯中的x和y的意义相反，因此推得公式的x，y要互换
                if (sp>=-(l-1)/2)&&(sp<=(l-1)/2)%若s’没有超出探测器的范围
                    isp = ceil(sp); sp0 = isp-sp;
                    proj_sp=(1-sp0)*R(isp+(l-1)/2+1,i+1)+sp0*R(isp+(l-1)/2,i+1);%根据s‘的值用他前后的射线投影为他插值，这里没用interp1，能够很大程度上加快运行速度
                    recon(x,y)=recon(x,y)+dbeta*(1/2*power(U(distance,(y-M2-1),(x-M2-1),r),2))*proj_sp;%重建图像上的像素点
                end
            end
        end
        str = ['重建进度：',num2str(round(i/359*100)),'%'];
        app.updateGUI(str);
        imshow(recon,[],'Parent',app.UIAxes3);
        drawnow
    end
end
