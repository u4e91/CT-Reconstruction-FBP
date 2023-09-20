%丁昊妍 520021910379
%% 解析法投影
imshow([],'Parent',app.UIAxes2);
imshow([],'Parent',app.UIAxes3);
imshow([],'Parent',app.UIAxes4);
xlabel([],'Parent',app.UIAxes2);
ylabel([],'Parent',app.UIAxes2);
xlabel([],'Parent',app.UIAxes4);
ylabel([],'Parent',app.UIAxes4);

M = 256;
randon_row=256;
radonImg_Computerized=double(zeros(randon_row,180));
E = app.LoadButton.UserData;
for i = 1:10 %计算每个椭圆
    rou = E(i,1);%灰度值
    a = E(i,2)*M/2;%椭圆横轴
    b = E(i,3)*M/2;%椭圆竖轴
    m = E(i,4)*M/2;%椭圆中心x
    n = E(i,5)*M/2;%椭圆中心y
    alpha = E(i,6)*pi/180;%椭圆横轴和x轴夹角
    for theta_degree =0:1:179 %计算每个角度
        theta = theta_degree*pi/180;
        a_theta_2 = a*a*((cos(theta-alpha))^2) + b*b*((sin(theta-alpha))^2);
        s = sqrt(m*m +n*n);
        
        if n == 0
            if m >= 0
                gamma =0;
            else
                gamma = pi;
            end
        else
            if m < 0
                gamma = atan(n/m) + pi;%二三象限
            else
                gamma = atan(n/m);%一四象限
            end
        end
        
        for j = 0:1:randon_row-1 %计算每个P
            t0 = j-0.5*randon_row;
            t = t0 -s*cos(gamma-theta);
            if t*t <= a_theta_2
                radonImg_Computerized(j+1,theta_degree+1) = radonImg_Computerized(j+1,theta_degree+1) + ((2*255*rou*a*b)/a_theta_2)*sqrt(a_theta_2 - t*t);
            end
        end
    end
    str = ['扫描进度：',num2str(round(i*10)),'%'];
    app.updateGUI(str);
    drawnow
end
radonImg_Computerized1 = imresize(radonImg_Computerized, [256,256]);
imshow(radonImg_Computerized1,[],'Parent',app.UIAxes2);
xlabel('\theta(degrees)','Parent',app.UIAxes2);
ylabel('x','Parent',app.UIAxes2);
drawnow
%% 傅里叶变换反投影
sensor_length = M;%投影器数量
filter_choose = app.filter_choose;
projection = radonImg_Computerized;
pixelsize(1) = M;
pixelsize(2) = M;
Img = fbp(projection,filter_choose,sensor_length,pixelsize, app);
imshow(Img,[],'Parent',app.UIAxes3);




