semilogy(snrdb,v1_0,'y-o','linewidth',2)
hold on
% semilogy(snrdb,v1_pi16,'r-o','linewidth',2)
% hold on
semilogy(snrdb,v1_pi8,'c-o','linewidth',2)
hold on
semilogy(snrdb,v1_pi4,'r-o','linewidth',2)
hold on
semilogy(snrdb,v2_0,'g-x','linewidth',2)
hold on
% semilogy(snrdb,v2_pi16,'m-x','linewidth',2)
% hold on
semilogy(snrdb,v2_pi8,'k-x','linewidth',2)
hold on
semilogy(snrdb,v2_pi4,'m-x','linewidth',2)
hold on
semilogy(snrdb,v0_0,'b-*','linewidth',2)

grid on
legend('Q2-0','Q2-pi8','Q2-pi4','Q1-0','Q1-pi8','Q1-pi4','1x1-0')
xlim([-4 30])
ylim([10^-8 10^0])
title('Q system BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label