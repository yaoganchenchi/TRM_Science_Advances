function data = loadMatData(filePath)
maskStruct = load(filePath);
data = getfield(maskStruct,char(fieldnames(maskStruct)));
clear maskStruct;
end