function cylinder = select_cylinders(cylinder,Ind)

Names = fieldnames(cylinder);
n = size(Names,1);
for i = 1:n
    cylinder.(Names{i}) = cylinder.(Names{i})(Ind,:);
end