function [R_gt,t_gt,data_x,data_y]=gen_data_anglenoise(num_inlier,num_outlier,noise_level)
num = num_inlier+num_outlier;

R_theta = (rand(1)*2-1)*1;
% R_gt = [cos(R_theta),sin(R_theta),0;-sin(R_theta),cos(R_theta),0;0,0,1];

v = rand(3,1)*2-1;
R_v = v./norm(v);
R_gt = rotationVectorToMatrix(R_theta*R_v);

scal = 100;
t_gt = (rand(3,1)*2-1)*scal;
T_gt = [R_gt,t_gt]';
fid = fopen('gt_syn.txt','w');
fprintf(fid,'%.10f %.10f %.10f %.10f\n',T_gt);

scale = 100;
data_x_inlier = (rand(3,num_inlier)*2-1)*scale;
data_y_inlier = R_gt*(data_x_inlier+t_gt);

data_x_outlier = (rand(3,num_outlier)*2-1)*scale;
data_y_outlier = (rand(3,num_outlier)*2-1)*scale;

data_x = [data_x_inlier,data_x_outlier]+normrnd(0,noise_level,3,num_inlier+num_outlier);
data_y = [data_y_inlier,data_y_outlier]+normrnd(0,noise_level,3,num_inlier+num_outlier);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save data
angle_noise1=0.1;
angle_noise2=0.3;
angle_noise3=0.5;

angle1=deg2rad(angle_noise1);
angle2=deg2rad(angle_noise2);
angle3=deg2rad(angle_noise3);

vv=rand(3,1)*2-1;
R_vv=vv./norm(vv);

R_noise1=rotationVectorToMatrix(angle1*R_vv);
R_noise2=rotationVectorToMatrix(angle2*R_vv);
R_noise3=rotationVectorToMatrix(angle3*R_vv);


v=data_x_inlier(:,1)-data_x_inlier(:,2);
v_p=v./norm(v);
v_q=R_gt*v_p;
v_q1=R_noise1*v_q;
v_q2=R_noise2*v_q;
v_q3=R_noise3*v_q;

v1=[v_p;v_q1]';
v2=[v_p;v_q2]';
v3=[v_p;v_q3]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('input_syn1.txt','w');
fprintf(fid,'%d\n',num);
fid=fopen('input_syn1.txt','a');
fprintf(fid,'%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n',v1);
for ii=1:num
    data_point=[data_y(:,ii);data_x(:,ii)]';
    fprintf(fid,'%.10f %.10f %.10f %.10f %.10f %.10f\n',data_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('input_syn2.txt','w');
fprintf(fid,'%d\n',num);
fid=fopen('input_syn2.txt','a');
fprintf(fid,'%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n',v2);
for ii=1:num
    data_point=[data_y(:,ii);data_x(:,ii)]';
    fprintf(fid,'%.10f %.10f %.10f %.10f %.10f %.10f\n',data_point);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('input_syn3.txt','w');
fprintf(fid,'%d\n',num);
fid=fopen('input_syn3.txt','a');
fprintf(fid,'%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n',v3);
for ii=1:num
    data_point=[data_y(:,ii);data_x(:,ii)]';
    fprintf(fid,'%.10f %.10f %.10f %.10f %.10f %.10f\n',data_point);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fid=fopen('input_syn4.txt','w');
% fprintf(fid,'%d\n',num);
% fid=fopen('input_syn4.txt','a');
% fprintf(fid,'%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n',v4);
% for ii=1:num
%     data_point=[data_y(:,ii);data_x(:,ii)]';
%     fprintf(fid,'%.10f %.10f %.10f %.10f %.10f %.10f\n',data_point);
% end

fclose('all');





end