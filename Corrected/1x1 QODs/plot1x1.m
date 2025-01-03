semilogy(snrdb,ber4c1,'m-*','linewidth',2)
hold on
semilogy(snrdb,ber4c2,'r-*','linewidth',2)
hold on
semilogy(snrdb,ber8c1,'b-o','linewidth',2)
hold on
semilogy(snrdb,ber8c2,'y-o','linewidth',2)
hold on
semilogy(snrdb,ber16c1,'c-x','linewidth',2)
hold on
semilogy(snrdb,ber16c2,'g-x','linewidth',2)

grid on
legend('PSK4,PSK4','PSK4,PSK4r','PSK8,PSK8','PSK8,PSK8r','PSK16,PSK16','PSK16,PSK16r')
xlim([-4 30])
ylim([10^-8 10^0])
title('1x1 system BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label