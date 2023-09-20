%丁昊妍 520021910379

imshow([],'Parent',app.UIAxes2);
imshow([],'Parent',app.UIAxes3);
imshow([],'Parent',app.UIAxes4);
xlabel([],'Parent',app.UIAxes2);
ylabel([],'Parent',app.UIAxes2);
xlabel([],'Parent',app.UIAxes4);
ylabel([],'Parent',app.UIAxes4);

P = app.LoadButton.UserData;

num_degrees = 180;
%确定发射器和接收器的位置
img_longest_edge = max(size(P,1), size(P,2));%图像最长边
sensor_length = ceil(img_longest_edge*sqrt(2))+2;%投影器数量
radonImg = zeros(sensor_length, num_degrees);
%projection_1d = zeros(sensor_length,1);


for theta_degree = 90:1:num_degrees+90-1
    theta = theta_degree*pi/180;
    projection_1d = zeros(sensor_length,1);
    %发射器中点
    mid_Trm_x = 0.5*(sensor_length - 2) * cos(theta + pi);
    mid_Trm_y = 0.5*(sensor_length - 2) * sin(theta + pi);
    %发射器
    Trm(1,1) = mid_Trm_x + 0.5*sensor_length*cos(0.5*pi+theta);
    Trm(2,1) = mid_Trm_y + 0.5*sensor_length*sin(0.5*pi+theta);
    for i = 2:1:sensor_length
        Trm(1,i) = Trm(1,i-1) - cos(0.5*pi+theta);
        Trm(2,i) = Trm(2,i-1) - sin(0.5*pi+theta);
    end

    %接收器中点
    mid_Rcv_x = 0.5*(sensor_length) * cos(theta);
    mid_Rcv_y = 0.5*(sensor_length) * sin(theta);
    %接收器
    Rcv(1,1) = mid_Rcv_x + 0.5*sensor_length*cos(0.5*pi+theta);
    Rcv(2,1) = mid_Rcv_y + 0.5*sensor_length*sin(0.5*pi+theta);
    for i = 2:1:sensor_length
        Rcv(1,i) = Rcv(1,i-1)-cos(0.5*pi+theta);
        Rcv(2,i) = Rcv(2,i-1)-sin(0.5*pi+theta);
    end
    
    %投影，这个角度下所有的投影值
    for i = 1:1:sensor_length
        projection_1d(i)=0;
        %计算相交比例alpha值范围
        [min_alpha, max_alpha, flag_intersection] = Caculate_range_alpha(i,Trm,Rcv,P);
        if flag_intersection == 0 %没有交点，或者仅有一个点的交点（而不是线段），所以投影为0
            continue;
        end
        %计算射线与图像相交的m,n范围
        [intersection_n_min, intersection_n_max, intersection_m_min, intersection_m_max] = Calculate_range_index_image(i, Trm, Rcv, P, min_alpha, max_alpha);

        %求解alpha_x, alpha_y
        if mod(theta_degree, 90)~=0
            %求alpha_x
            sizeof_alpha_x = round(intersection_n_max-intersection_n_min+1);
            alpha_x = zeros(sizeof_alpha_x,1);
            for k = 0:1:sizeof_alpha_x-1
                x_position = intersection_n_min + k - size(P,2)/2;
                alpha_x(k+1) = (x_position -Trm(1,i))/(Rcv(1,i)-Trm(1,i));
            end
            %求alpha_y
            sizeof_alpha_y = round(intersection_m_max-intersection_m_min+1);
            alpha_y = zeros(sizeof_alpha_y,1);
            for k = 0:1:sizeof_alpha_y-1
                y_position = size(P,1) / 2 - (intersection_m_min + k);
                alpha_y(k+1) = (y_position -Trm(2,i))/(Rcv(2,i)-Trm(2,i));
            end

            %合并
            sizeof_alpha = sizeof_alpha_x + sizeof_alpha_y;
            alpha = [alpha_x; alpha_y];
            alpha = sort(alpha);%, alpha + sizeof_alpha);%升序排序
        elseif mod(theta_degree, 180) == 0%平行x轴
            %求alpha_x
            sizeof_alpha_x = round(intersection_n_max-intersection_n_min+1);
            alpha_x = zeros(sizeof_alpha_x,1);
            for k = 0:1:sizeof_alpha_x-1
                x_position = intersection_n_min + k - size(P,2)/2;
                alpha_x(k+1) = (x_position -Trm(1,i))/(Rcv(1,i)-Trm(1,i));
            end
            %求alpha
            sizeof_alpha = sizeof_alpha_x;
            alpha = alpha_x;
            alpha = sort(alpha);%, alpha + sizeof_alpha);
        else
            %求alpha_y
            sizeof_alpha_y = round(intersection_m_max-intersection_m_min+1);
            alpha_y = zeros(sizeof_alpha_y,1);
            for k = 0:1:sizeof_alpha_y-1
                y_position = size(P,1) / 2 - (intersection_m_min + k);
                alpha_y(k+1) = (y_position -Trm(2,i))/(Rcv(2,i)-Trm(2,i));
            end
            sizeof_alpha = sizeof_alpha_y;
            alpha = alpha_y;
            alpha = sort(alpha);%, alpha + sizeof_alpha);
        end
        if i == 61
            grey = 0;
        end
        projection_1d(i) = Calculate_projection_oneValue(i, Trm, Rcv, P, sizeof_alpha, alpha);
    end

    for m = 1:1:size(radonImg,1)
        radonImg (m, theta_degree-90+1) = projection_1d(m);
    end
    str = ['扫描进度：',num2str(round((theta_degree-90+1)/180*100)),'%'];
    app.updateGUI(str);
    drawnow
    
    radonImg1 = imresize(radonImg, [256,256]);
    imshow(radonImg1,[],'Parent',app.UIAxes2);
    xlabel('\theta(degrees)','Parent',app.UIAxes2);
    ylabel('x','Parent',app.UIAxes2);
end
radonImg1 = imresize(radonImg, [256,256]);
imshow(radonImg1,[],'Parent',app.UIAxes2);
xlabel('\theta(degrees)','Parent',app.UIAxes2);
ylabel('x','Parent',app.UIAxes2);
drawnow
%% 傅里叶变换反投影
filter_choose = app.filter_choose;
projection = radonImg;
pixelsize(1) = size(P,1);
pixelsize(2) = size(P,2);
Img = fbp(projection,filter_choose,sensor_length,pixelsize, app);
imshow(Img,[],'Parent',app.UIAxes3);




function [min_alpha, max_alpha, flag_intersection]  = Caculate_range_alpha(i,Trm,Rcv,P)%计算alpha范围
    %确定图像位置
    P_col = size(P,2);
    P_row = size(P,1);
    img_x_min = -0.5*P_col;
    img_y_min = -0.5*P_row;
    img_x_max = img_x_min + P_col;
    img_y_max = img_x_min + P_row;
    min_alpha_x = 1.0;
    min_alpha_y = 1.0;
    max_alpha_x = 0.0;
    max_alpha_y = 0.0;
    if (Rcv(1,i)-Trm(1,i))~=0 %不竖直
        min_alpha_x = min((img_x_min - Trm(1,i)) / (Rcv(1,i) - Trm(1,i)),...
            (img_x_max - Trm(1,i)) / (Rcv(1,i) - Trm(1,i)));
        max_alpha_x = max((img_x_min - Trm(1,i)) / (Rcv(1,i) - Trm(1,i)),...
            (img_x_max - Trm(1,i)) / (Rcv(1,i) - Trm(1,i)));
    end
    if ((Rcv(2,i) - Trm(2,i)) ~= 0)%不平行
        min_alpha_y = min((img_y_min - Trm(2,i)) / (Rcv(2,i) - Trm(2,i)),...
            (img_y_max - Trm(2,i)) / (Rcv(2,i) - Trm(2,i)));
        max_alpha_y = max((img_y_min - Trm(2,i)) / (Rcv(2,i) - Trm(2,i)),...
            (img_y_max - Trm(2,i)) / (Rcv(2,i) - Trm(2,i)));
    end

    min_alpha = max([min_alpha_x, min_alpha_y, 0.0]);
    max_alpha = min([max_alpha_x, max_alpha_y, 1.0]);

    flag_intersection = 1;
    if(min_alpha >= max_alpha)
        flag_intersection = 0;
    end
end


function [intersection_n_min, intersection_n_max, intersection_m_min, intersection_m_max] = Calculate_range_index_image(i, Trm, Rcv, P, min_alpha, max_alpha)
    %确定图像位置
    P_col = size(P,2);
    P_row = size(P,1);
    img_x_min = -0.5*P_col;
    img_y_min = -0.5*P_row;
    img_x_max = img_x_min + P_col;
    img_y_max = img_x_min + P_row;

    if Trm(1,i) < Rcv(1,i)
        intersection_n_min = (P_col) - (img_x_max - min_alpha * (Rcv(1,i) - Trm(1,i)) - Trm(1,i));
        intersection_n_max = -1 + (Trm(1,i) + max_alpha * (Rcv(1,i) - Trm(1,i)) - img_x_min);
    else
        intersection_n_min = (P_col) - (img_x_max - max_alpha * (Rcv(1,i) - Trm(1,i)) - Trm(1,i));
        intersection_n_max = -1 + (Trm(1,i) + min_alpha * (Rcv(1,i) - Trm(1,i)) - img_x_min);
    end


    if Trm(2,i) < Rcv(2,i)
        intersection_m_min = (P_col) - (Trm(2,i) + max_alpha * (Rcv(2,i) - Trm(2,i)) - img_y_min);
        intersection_m_max = img_y_max - min_alpha * (Rcv(2,i) - Trm(2,i)) - Trm(2,i) - 1;
    else
        intersection_m_min = (P_col) - (Trm(2,i) + min_alpha * (Rcv(2,i) - Trm(2,i)) - img_y_min);
        intersection_m_max = img_y_max - max_alpha * (Rcv(2,i) - Trm(2,i)) - Trm(2,i) - 1;
    end


end

function projection = Calculate_projection_oneValue(i, Trm, Rcv, P, sizeof_alpha, alpha)
    total_TR_length = sqrt((Trm(1,i)-Rcv(1,i))^2+(Trm(2,i)-Rcv(2,i))^2);
    projection = 0.0;
    P_col = size(P,2);
    P_row = size(P,1);
    for k=2:1:sizeof_alpha
        length_little_route = (alpha(k)-alpha(k-1))*total_TR_length;
        mid_alpha = 0.5*(alpha(k)+alpha(k-1));
        M_pixel(1) = Trm(1,i) + mid_alpha*(Rcv(1,i)-Trm(1,i));
        M_pixel(2) = Trm(2,i) + mid_alpha*(Rcv(2,i)-Trm(2,i));

        n_pixel = floor(M_pixel(1) + 0.5 * P_col);
        n_pixel = max(0, n_pixel);
        %n_pixel = min(n_pixel, P_row - 1);
        n_pixel = min(n_pixel, P_col - 1);

        m_pixel = floor(0.5 * P_row - M_pixel(2));
        m_pixel = max(0, m_pixel);
        %m_pixel = min(m_pixel, P_col - 1);
        m_pixel = min(m_pixel, P_row - 1);

        grayscale = P(m_pixel+1, n_pixel+1);
        projection = projection + length_little_route * grayscale;
    end
end