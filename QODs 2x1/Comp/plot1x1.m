semilogy(snrdb,QOD1x1_0,'m-*','linewidth',2)
hold on
semilogy(snrdb,QOD1x1_pi4,'r-*','linewidth',2)
hold on
semilogy(snrdb,QOD2x1_00,'b-x','linewidth',2)
hold on
semilogy(snrdb,QOD2x1_pi16,'y-x','linewidth',2)
hold on
semilogy(snrdb,QOD2x1_pi8,'c-x','linewidth',2)
hold on
semilogy(snrdb,QOD2x1_pi4,'g-x','linewidth',2)
hold on
semilogy(snrdb,QOD2x1_Sympi4,'k-o','linewidth',2)

grid on
legend('QOD1x1-0','QOD1x1-pi4','QOD2x1-0','QOD2x1-pi16','QOD2x1-pi8','QOD2x1-pi4','QOD2x1-Sympi4')
xlim([-4 30])
ylim([10^-8 10^0])
title('QPSK comparion BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label