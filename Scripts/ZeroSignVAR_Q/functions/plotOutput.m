function plotOutput( y, identifiedModels, output, nVars, t, nModelsFound, opt, ctmModel )

    if opt.isImpulsePlots
        % Plot IRFs in single plots
        xAxis = 0:opt.nImpulseHorizon;
        xFill = [xAxis,fliplr(xAxis)];
        yInnerFill = [output.irf16',fliplr(output.irf84')];
        yOuterFill = [output.irf05',fliplr(output.irf95')];
        if opt.nPlotRandomModels
            randomModels = randi(size(identifiedModels,2),opt.nPlotRandomModels,1);
        end
        jj = 1;
        for ii = 1:nVars^2
            if mod(ii,nVars) == 0
                shock = opt.lShocks{nVars};
            else
                shock = opt.lShocks{mod(ii,nVars)};
            end
            if ~isempty(shock)
                fig = figure(ii);
                hold on
                    fill(xFill,yOuterFill(ii,:),[0.8 0.8 0.8],'edgecolor','none');
                    fill(xFill,yInnerFill(ii,:),[0.75 0.75 0.75],'edgecolor','none');
                    if opt.isPlotCTM
                        plot(0:opt.nImpulseHorizon,output.irfCTM(1:opt.nImpulseHorizon+1,ii),'k:');
                    end
                    plot(0:opt.nImpulseHorizon,output.irfMedian(1:opt.nImpulseHorizon+1,ii),'k-');
                    plot(0:opt.nImpulseHorizon,zeros(opt.nImpulseHorizon+1)','k-');
                    if opt.nPlotRandomModels    
                        for kk = 1:opt.nPlotRandomModels
                            plot(0:opt.nImpulseHorizon,identifiedModels(randomModels(kk)).irf(1:opt.nImpulseHorizon+1,ii),'--*r');
                        end 
                    end
                    xlabel('Horizon');
                    xlim([0 opt.nImpulseHorizon])
                hold off
                if opt.isPlotTitle
                    title([shock ': ' opt.lVars{jj}]);
                end                
                saveas(fig, strcat(opt.modelPath,opt.modelName,'/imgs/',shock,'-',opt.lVars{jj}),'png');
            end
            if mod(ii,nVars) == 0
                jj = jj+1;
            end
        end
        close all;
    end
     
    if opt.isIRFtable == 1
        save(strcat(opt.modelPath,opt.modelName,'/tables/irfs.mat'), '-struct', 'output', 'irfCTM', 'irfMedian', 'irf95', 'irf84', 'irf16', 'irf05');
    end
    
    
    if opt.isStructShockPlots

        startDate = datenum(opt.startDate);
        endDate = datenum(opt.endDate);
        xData = linspace(startDate,endDate,t-opt.nLags);
        
%         xAxis = 1:550;
%         xFill = [xAxis,fliplr(xAxis)];
%         yInnerFill = [output.structShock16',fliplr(output.structShock84')];
%         yOuterFill = [output.structShock05',fliplr(output.structShock95')];        
        for ii = 1:nVars
            fig = figure(ii);
            hold on
%                 fill(xFill,yOuterFill(ii,:),[0.8 0.8 0.8],'edgecolor','none');
%                 fill(xFill,yInnerFill(ii,:),[0.75 0.75 0.75],'edgecolor','none');
                plot(xData,output.structShockMedian(:,ii),'k-');
                datetick('x','yyyy');
            hold off
            if opt.isPlotTitle
                title(opt.lShocks{ii});
            end
            saveas(fig, strcat(opt.modelPath,opt.modelName,'/imgs/shocks-',opt.lShocks{ii}),'png');
        end
        close all;
    end
    
    if opt.isStructShockTable
        output.structShockCTM = ctmModel.structShocks;
        save(strcat(opt.modelPath,opt.modelName,'/tables/StructShocks.mat'), '-struct', 'output', 'structShockCTM', 'structShockMedian', 'structShock95', 'structShock84', 'structShock16', 'structShock05');
    end    
    
    if opt.isFevdTexTable == 1
        for ii = 1:nVars
            fid = fopen(strcat(opt.modelPath,opt.modelName,'/tables/fevd-',opt.lVars{ii},'.tex'), 'w+');
            fprintf(fid, '\\documentclass{article}\n');
            fprintf(fid, '\\usepackage{multirow}\n');
            fprintf(fid, '\\begin{document}\n');
            fprintf(fid, '\\begin{table}\n');
            fprintf(fid, '\t\\footnotesize\n');
            fprintf(fid, ['\t\\caption{FEVD of ' opt.lVars{ii} ' for ' opt.modelName '}\n']);
            fprintf(fid, '\t%s%s%s\n', '\begin{tabular}{', repmat('c',1,nVars+1), '}');
            fprintf(fid, '\t\tHorizon');
            fprintf(fid, ' & %s', opt.lShocks{1:nVars});
            fprintf(fid, '\\\\ \\hline\n');
            for jj = 1:opt.nImpulseHorizon
                fprintf(fid, '\t\t\\multirow{2}{*}{%d}', jj);
                fprintf(fid, ' & %2.2f', output.fevdMedian(jj,(ii-1)*nVars+1:ii*nVars));
                fprintf(fid, '\\\\\n');
                fprintf(fid, '\t\t');
                temp = [output.fevd16(jj,(ii-1)*nVars+1:ii*nVars); output.fevd84(jj,(ii-1)*nVars+1:ii*nVars)];
                fprintf(fid, ' & (%2.2f, %2.2f)', temp(:)');
                clear temp;
                fprintf(fid, '\\\\\n');
            end
            fprintf(fid, '\t\\end{tabular}\n');
            fprintf(fid, strcat('\\label{tab:fevd-', opt.modelName, '-', opt.lVars{ii}, '}\n'));
            fprintf(fid, '\\end{table}');
            fprintf(fid, '\\end{document}');
            fclose(fid);        
        end       
    end
    
    if opt.isFevdTable == 1
        save(strcat(opt.modelPath,opt.modelName,'/tables/fevd.mat'), '-struct', 'output', 'fevdCTM', 'fevdMedian', 'fevd95', 'fevd84', 'fevd16', 'fevd05');
    end
    
    if opt.isFevdPlots==1 && opt.isPlotCTM==0
        for ii = 1:nVars
            fig = figure(ii);
            area(output.fevdCTM(:,(ii-1)*nVars+1:ii*nVars));
            colormap gray;
            set(gca,'YLim',[0,100]);
            legend(opt.lShocks);
            xlabel('Horizon');
            if opt.isPlotTitle
                title(opt.lVars{ii});
            end
            saveas(fig, strcat(opt.modelPath,opt.modelName,'/imgs/fevd-',opt.lVars{ii}),'png');
        end
        close all;
    end
      
    if opt.isFevdPlots==1 && opt.isPlotCTM==1
        for ii = 1:nVars
            fig = figure(ii);
            area(output.fevdCTM(:,(ii-1)*nVars+1:ii*nVars));
            colormap gray;
            set(gca,'YLim',[0,100]);
            legend(opt.lShocks);
            xlabel('Horizon');
            if opt.isPlotTitle
                title(opt.lVars{ii});
            end
            saveas(fig, strcat(opt.modelPath,opt.modelName,'/imgs/fevdCTM-',opt.lVars{ii}),'png');
        end
        close all;
        for ii = 1:nVars
            fig = figure(ii);
            area(output.fevdMedian(:,(ii-1)*nVars+1:ii*nVars));
            colormap gray;
            set(gca,'YLim',[0,100]);
            legend(opt.lShocks);
            xlabel('Horizon');
            if opt.isPlotTitle
                title(opt.lVars{ii});
            end
            saveas(fig, strcat(opt.modelPath,opt.modelName,'/imgs/fevdPWM-',opt.lVars{ii}),'png');
        end
        close all;
    end
    
    
    if opt.isHistDecompPlots

        startDate = datenum(opt.startDate);
        endDate = datenum(opt.endDate);
        xData = linspace(startDate,endDate,t-opt.nLags);
        
        for ii = 1:nVars
            
            tempMatrixPos = output.HistDecompMedian(:,(ii-1)*nVars+1:ii*nVars);
            tempMatrixNeg = output.HistDecompMedian(:,(ii-1)*nVars+1:ii*nVars);
            tempMatrixPos(tempMatrixPos<0) = NaN;
            tempMatrixNeg(tempMatrixNeg>0) = NaN;
            
            fig = figure(ii);
            hold on
                tempPlot1 = bar(xData, tempMatrixPos, .9, 'stack', 'EdgeColor','none');
                bar(xData, tempMatrixNeg, .9, 'stack', 'EdgeColor','none');
                %plot(xData,y(opt.nLags+1:end,ii),'k-');
                tempPlot2 = plot(xData,y(identifiedModels(1).nLags+1:end,ii)-output.BaseProjectMedian(:,ii),'k-', 'LineWidth', 1);
                datetick('x','yyyy');
            hold off
            if opt.isPlotTitle
                title(opt.lVars{ii});
            end
            tempLeg = [opt.lShocks, 'Deviation'];
            legend([tempPlot1 tempPlot2], tempLeg, 'Location', 'NorthEastOutside')
            saveas(fig, strcat(opt.modelPath,opt.modelName,'/imgs/HistDecompPWM-',opt.lVars{ii}),'pdf');
        end
        close all;
        
        if opt.isPlotCTM==1
        for ii = 1:nVars
            
            tempMatrixPos = ctmModel.HistDecomp(:,(ii-1)*nVars+1:ii*nVars);
            tempMatrixNeg = ctmModel.HistDecomp(:,(ii-1)*nVars+1:ii*nVars);
            tempMatrixPos(tempMatrixPos<0) = NaN;
            tempMatrixNeg(tempMatrixNeg>0) = NaN;
            
            fig = figure(ii);
            hold on
                tempPlot1 = bar(xData, tempMatrixPos, .9, 'stack', 'EdgeColor','none');
                bar(xData, tempMatrixNeg, .9, 'stack', 'EdgeColor','none');
                %plot(xData,y(opt.nLags+1:end,ii),'k-');
                tempPlot2 = plot(xData,y(identifiedModels(1).nLags+1:end,ii)-ctmModel.BaseProject(:,ii),'k-', 'LineWidth', 1);
                datetick('x','yyyy');
            hold off
            if opt.isPlotTitle
                title(opt.lVars{ii});
            end
            tempLeg = [opt.lShocks, 'Deviation'];
            legend([tempPlot1 tempPlot2], tempLeg, 'Location', 'NorthEastOutside')
            saveas(fig, strcat(opt.modelPath,opt.modelName,'/imgs/HistDecompCTM-',opt.lVars{ii}),'pdf');
        end
        close all;
        end
    end
    
    if opt.isHistDecompTable
        save(strcat(opt.modelPath,opt.modelName,'/tables/HistDecompPWM.mat'), '-struct', 'output', 'HistDecompMedian', 'HistDecomp95', 'HistDecomp84', 'HistDecomp16', 'HistDecomp05', 'BaseProjectMedian');
        save(strcat(opt.modelPath,opt.modelName,'/tables/HistDecompCTM.mat'), '-struct', 'ctmModel', 'HistDecomp', 'BaseProject');
    end    
end