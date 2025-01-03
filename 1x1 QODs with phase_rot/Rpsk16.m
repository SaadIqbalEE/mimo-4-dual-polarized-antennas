clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The PSK16 modulation %%%  

tic;        %Timer Started

%Two modulation Schemes QPSK & QPSKr
% PSK16 <--- Z1
% PSK16r<--- Z2

l= floor(100000/8)*8;            %No. of Bits @Note: should be divisble by nbits

snrdb=-4:2:30;      %Array of SNR
snrx= 10.^(-snrdb/20);
ber16c1 = zeros(1,length(snrdb));
ber16c2 = zeros(1,length(snrdb));


bits=round(rand(1,l));  %Bits generation
bit = 0;                %Resultant bits

zpsk16 = exp((1i*(pi/8)*[0:15])+(0));    %psk16 "Unit Power"
zpsk16r  = exp((1i*(pi/8)*[0:15])+(pi/16));    %psk16 rotated "Unit Power"

M1 = length(zpsk16);
M2 = length(zpsk16r);

nbits1 = log2(M1);
sym1=zeros(1,l/nbits1);   %Symbol sent -> QPSK

nbits2 = log2(M2);
sym2=zeros(1,l/nbits2);   %Symbol sent -> QPSK

for len = nbits1:nbits1:l
    switch (bi2de([bits(len) bits(len-1) bits(len-2) bits(len-3)]))
        case 0
            sym1(len/nbits1) = zpsk16(1);
        case 1
            sym1(len/nbits1) = zpsk16(2);         
        case 2
            sym1(len/nbits1) = zpsk16(3);            
        case 3
             sym1(len/nbits1) = zpsk16(4);
        case 4
            sym1(len/nbits1) = zpsk16(5);
        case 5
            sym1(len/nbits1) = zpsk16(6);         
        case 6
            sym1(len/nbits1) = zpsk16(7);            
        case 7
             sym1(len/nbits1) = zpsk16(8);
        case 8
            sym1(len/nbits1) = zpsk16(9);
        case 9
            sym1(len/nbits1) = zpsk16(10);         
        case 10
            sym1(len/nbits1) = zpsk16(11);            
        case 11
             sym1(len/nbits1) = zpsk16(12);
        case 12
            sym1(len/nbits1) = zpsk16(13);
        case 13
            sym1(len/nbits1) = zpsk16(14);         
        case 14
            sym1(len/nbits1) = zpsk16(15);            
        case 15
             sym1(len/nbits1) = zpsk16(16);
    end
end

for len = nbits2:nbits2:l
    switch (bi2de([bits(len) bits(len-1) bits(len-2) bits(len-3)]))
        case 0
            sym2(len/nbits1) = zpsk16r(1);
        case 1
            sym2(len/nbits1) = zpsk16r(2);         
        case 2
            sym2(len/nbits1) = zpsk16r(3);            
        case 3
            sym2(len/nbits1) = zpsk16r(4);
        case 4
            sym2(len/nbits1) = zpsk16r(5);
        case 5
            sym2(len/nbits1) = zpsk16r(6);         
        case 6
            sym2(len/nbits1) = zpsk16r(7);            
        case 7
            sym2(len/nbits1) = zpsk16r(8);
        case 8
            sym2(len/nbits1) = zpsk16r(9);
        case 9
            sym2(len/nbits1) = zpsk16r(10);         
        case 10
            sym2(len/nbits1) = zpsk16r(11);            
        case 11
            sym2(len/nbits1) = zpsk16r(12);
        case 12
            sym2(len/nbits1) = zpsk16r(13);
        case 13
            sym2(len/nbits1) = zpsk16r(14);         
        case 14
            sym2(len/nbits1) = zpsk16r(15);            
        case 15
            sym2(len/nbits1) = zpsk16r(16);
    end
end

out = zeros(1,l/nbits1);     %Received symbols array

SD = 1/sqrt(2);
loop_Cond = sym1;

for Case=1:2
    for n=1:length(snrdb)
        for loop=1:2:size(sym1,2)        %Two symbols at every iteration.

            Hx11=normrnd(0,SD,[1 1]);
            Hy11=normrnd(0,SD,[1 1]);
            H1=(Hx11+1i*Hy11);

            Hx22=normrnd(0,SD,[1 1]);
            Hy22=normrnd(0,SD,[1 1]);
            H2=(Hx22+1i*Hy22);

            w1_0 = snrx(n)*(1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));        %Noise Profile RH
            w1_1 = snrx(n)*(1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));        %Noise Profile RV

            if Case==1
                q1=sym1(loop);                   % First symbol
                q2=sym1(loop+1);                 % Second symbol
            elseif Case==2
                q1=sym1(loop);                   % First symbol
                q2=sym2(loop+1);                 % Second symbol
            end

            RH = q1*H1 - q2*conj(H2) + w1_0;      %Received at RH            
            RV = q1*H2 + q2*conj(H1) + w1_1;      %Received at RV

    %%%%%%%% Received Signal at Decoder %%%%%%%%%%%

            if Case==1
                lenz1 = size(zpsk16);
                lenz2 = size(zpsk16);
            elseif Case==2
                lenz1 = size(zpsk16);
                lenz2 = size(zpsk16r);
            end
            
            const=-100000;
            for z1=1:lenz1(2)
                for z2=1:lenz2(2)

                    C1 = zpsk16(z1);
                    
                    if Case==1  
                        C2 = zpsk16(z2);
                    elseif Case==2
                        C2 = zpsk16r(z2);
                    end

                    DR = (real(conj(RH)*(C1*H1-C2*conj(H2)) + RV*(conj(C1)*conj(H2)+conj(C2)*H1)));       %Decoding Symbol
                    if (DR > const)
                        temp = [z1 z2];
                        const= DR;
                    end
                end
            end
            
            if Case==1  
                out(loop) = zpsk16(temp(1));
                out(loop+1) = zpsk16(temp(2));
            elseif Case==2
                out(loop) = zpsk16(temp(1));
                out(loop+1) = zpsk16r(temp(2));
            end            
        end
        
        if Case==1  
            for k=1:1:l/nbits1
                switch out(k)
                    case zpsk16(1)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0 0 0];
                    case zpsk16(2)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0 0 1];                
                    case zpsk16(3)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0 1 0];                
                    case zpsk16(4)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0 1 1];
                    case zpsk16(5)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1 0 0];
                    case zpsk16(6)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1 0 1];                
                    case zpsk16(7)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1 1 0];                
                    case zpsk16(8)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1 1 1];
                    case zpsk16(9)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0 0 0];
                    case zpsk16(10)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0 0 1];                
                    case zpsk16(11)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0 1 0];                
                    case zpsk16(12)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0 1 1];
                    case zpsk16(13)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1 0 0];
                    case zpsk16(14)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1 0 1];                
                    case zpsk16(15)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1 1 0];                
                    case zpsk16(16)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1 1 1];
                end
            end
            ber16c1(n)=(l-sum(bits==bit))/l;
        elseif Case==2
            for k=1:1:l/nbits1
                if (rem(k,2) == 1)
                    switch out(k)
                        case zpsk16(1)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0 0 0];
                        case zpsk16(2)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0 0 1];                
                        case zpsk16(3)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0 1 0];                
                        case zpsk16(4)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0 1 1];
                        case zpsk16(5)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1 0 0];
                        case zpsk16(6)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1 0 1];                
                        case zpsk16(7)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1 1 0];                
                        case zpsk16(8)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1 1 1];
                        case zpsk16(9)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0 0 0];
                        case zpsk16(10)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0 0 1];                
                        case zpsk16(11)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0 1 0];                
                        case zpsk16(12)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0 1 1];
                        case zpsk16(13)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1 0 0];
                        case zpsk16(14)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1 0 1];                
                        case zpsk16(15)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1 1 0];                
                        case zpsk16(16)
                           bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1 1 1];
                    end
                else
                    switch out(k)
                        case zpsk16r(1)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 0 0 0];
                        case zpsk16r(2)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 0 0 1];                
                        case zpsk16r(3)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 0 1 0];                
                        case zpsk16r(4)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 0 1 1];
                        case zpsk16r(5)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 1 0 0];
                        case zpsk16r(6)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 1 0 1];                
                        case zpsk16r(7)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 1 1 0];                
                        case zpsk16r(8)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 1 1 1];
                       case zpsk16r(9)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 0 0 0];
                        case zpsk16r(10)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 0 0 1];                
                        case zpsk16r(11)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 0 1 0];                
                        case zpsk16r(12)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 0 1 1];
                        case zpsk16r(13)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 1 0 0];
                        case zpsk16r(14)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 1 0 1];                
                        case zpsk16r(15)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 1 1 0];                
                        case zpsk16r(16)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 1 1 1];
                    end
                end
            end
            ber16c2(n)=(l-sum(bits==bit))/l;
        end  
    end
end
    
semilogy(snrdb,ber16c1,'m-*','linewidth',2)
hold on
semilogy(snrdb,ber16c2,'r-*','linewidth',2)
grid on
legend('PSK8,PSK8','PSK8,PSK8r')
xlim([-4 30])
ylim([10^-8 10^0])
title('1x1 system BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label
toc 


        