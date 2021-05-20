function saveModelInfos( identifiedModels, S, Z, nModelsFound, nModelsChecked, nVars, opt )
    fileID = fopen(strcat(opt.modelPath,opt.modelName,'/infos.txt'),'w');
    fprintf(fileID,'%40s %25s\n','Model name:',opt.modelName);
	fprintf(fileID,'%40s %25s\n','Estimation Method:',opt.estimationMethod);
    fprintf(fileID,'%40s %25s\n','number of lags:',int2str(opt.nLags));
    if ~strcmp(opt.estimationMethod,'OLS')
        fprintf(fileID,'%40s %25s\n','draws to generate posterior:',int2str(opt.nDrawsFromBvar));
        fprintf(fileID,'%40s %25s\n','(max) number of draws from posterior:',int2str(opt.nModelDraws));
        fprintf(fileID,'%40s %25s\n','number of transformations per draw:',int2str(opt.nTransformationsPerDraw));        
    end
    if opt.keepAllValid
        fprintf(fileID,'%40s %25s\n','Transform only until match is found:','NO');
    else
        fprintf(fileID,'%40s %25s\n','Transform only until match is found:', 'YES');
    end
    fprintf(fileID,'%40s %25s\n','overall number of models checked:',int2str(nModelsChecked));
    fprintf(fileID,'%40s %25s\n','Models found:',int2str(nModelsFound));

    fclose(fileID);
end