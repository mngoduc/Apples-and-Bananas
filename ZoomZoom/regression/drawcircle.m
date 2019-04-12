function h = drawcircle(x,y,r,color,fill)
step=0.01;
ang=0:step:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
x_data = x+xp;
y_data = y+yp;
if fill == 1
rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1,1],'FaceColor',color,'LineStyle','none');
else
h = plot(x_data,y_data,color,'LineWidth',1.1);
end
hold on;
end

