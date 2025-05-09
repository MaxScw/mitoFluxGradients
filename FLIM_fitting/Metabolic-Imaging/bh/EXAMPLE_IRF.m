% execute this example by hitting F5 in Matlab
%% read an IRF
fn= 'IRF.sdt';
fp= fullfile( fn);

sp = bh_readsetup(fp);
disp(sp);   % show header
disp('no of data blocks in the file');
disp(sp.no_of_data_blocks);
blk_info = bh_blockinfo(sp);
for ii=1:sp.no_of_data_blocks
    dat(ii).blk_info= blk_info{ii};
    disp(['#' num2str(ii) ': ' dat(ii).blk_info.mode]);
    dat(ii).measdesc= bh_getmeasdesc(sp,ii);    
    dat(ii).data= bh_getdatablock(sp, ii);  % read data block 
    % calc. time resolution:
    dat(ii).dt= double(dat(ii).measdesc.tac_r)/double(dat(ii).measdesc.tac_g)/double(dat(ii).measdesc.adc_re);
    dat(ii).time_ns= double( [0:(dat(ii).measdesc.adc_re-1)]')* dat(ii).dt*1e9 ; % time axis in ns
end
%% Plot
figure; 
for ii=1:sp.no_of_data_blocks
    semilogy( dat(ii).time_ns, dat(ii).data); hold on;
end
grid on;
xlabel('time / ns'); ylabel('counts'); title('www.becker-hickl.com');