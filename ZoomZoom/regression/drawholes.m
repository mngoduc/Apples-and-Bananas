function drawholes(hole_data,numholes,color,fill)
for i=1:numholes
    drawcircle(hole_data(i,1),hole_data(i,2),hole_data(i,3),color,fill);
end



