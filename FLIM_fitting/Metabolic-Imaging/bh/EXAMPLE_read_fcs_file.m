% execute this example by hitting F5 in Matlab
%% Example read FCS data
fn= 'fcs.sdt';
fp= fullfile( fn);

sp = bh_readsetup(fp);
disp('no of data blocks in the file');
disp(sp.no_of_data_blocks);
blk_info = bh_blockinfo(sp)
for ii=1:sp.no_of_data_blocks
    dat(ii).blk_info= blk_info{ii};
    disp(['#' num2str(ii) ': ' dat(ii).blk_info.mode]);
    % assign correct measurement description block:
    ii_mb= dat(ii).blk_info.meas_desc_block_no + 1 ;
    dat(ii).measdesc= bh_getmeasdesc(sp,ii_mb);
    dat(ii).data= bh_getdatablock(sp, ii);
end
%% Plot
figure;
subplot(1,2,1);
single_ind= [1 4];
for ii=1:numel(single_ind)
    jj= single_ind(ii); 
    % calc. time resolution:
    dat(jj).dt= double(dat(jj).measdesc.tac_r)/double(dat(jj).measdesc.tac_g)/double(dat(jj).measdesc.adc_re);
    dat(jj).time_ns= double( [0:(dat(jj).measdesc.adc_re-1)]')* dat(jj).dt*1e9 ; % time axis in ns
    
    ph= semilogy( dat( jj ).time_ns, dat( jj ).data);
    set(ph, 'DisplayName', ['Decay, M' num2str(dat(jj).measdesc.fcs_mod+1) ] );
    hold all
end
grid on; axis tight; legend show
xlabel('time / ns'); ylabel('counts'); title('www.becker-hickl.com');

subplot(1,2,2);
fcs_ii=[2 3 5] 
for ii=1:numel(fcs_ii)
    jj= fcs_ii(ii); 
    ph= semilogx( dat( jj ).data(:,1), dat( jj ).data(:,2));
    set(ph, 'DisplayName', ['F(C)CS M' num2str(dat(jj).measdesc.fcs_mod+1) ... 
        ' vs. M' num2str( dat(jj).measdesc.fcs_cross_mod+1) ] );
    hold all
end
legend show
xlim([1e-2 1e5]);
ylim([.99 1.1]); grid on; title('www.becker-hickl.com');
xlabel('time / \mus'); ylabel('correlation factor'); 
