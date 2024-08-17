clear;
clc;
% Time-domain simulation to extract PFR power
runtime = 120 ;
simulationData = zeros(1+runtime/0.02);
load_system('D:\Desktop\SimulinkModel\QFproof.slx'); % Loading
num_iterations= 500 ;
result = zeros(10, num_iterations);
for ma = 1:1:num_iterations
    H = 3 + (8 - 2) * rand;
    k = 18 + (25 - 18) * rand;
    F= 0.15 + (0.4 - 0.15) * rand;
    F = 0;
    P= -(0.05 + (0.15 - 0.05) * rand);
    
    T= 6 + (14 - 6) * rand;
    result(1,ma) = H;
    result(2,ma) = k;
    result(3,ma) = F;
    result(4,ma) = T;
    result(5,ma) = P;
 
    set_param('QFproof/Inertia', 'Denominator', ['[', num2str(2*H), ',1]']); 
    set_param('QFproof/G', 'Numerator', ['[', num2str(F*T), ',1]']) 
    set_param('QFproof/G','Denominator', ['[', num2str(T), ',1]']) 
    set_param('QFproof/Gain','Gain', [num2str(-k)]) 
    set_param('QFproof/deltaP','After', [num2str(P)]) 

    % Run
    simout = sim('QFproof.slx',[0 runtime]);
    Frequncey = simout.get('Frequency');
    Power = simout.get('Power');

    [min_value, min_index] = min(Frequncey);
    tnadir = Power(min_index,1);
    P_gen = Power(min_index,2);

    fmax = min(Frequncey);
    fmax = abs(fmax(2));
    
    % QF
    x = abs(2 * fmax / (P / (2 * H )));
    P_QF = 2 * T * fmax * k * (((x + 2 * T) / (2 * T * x)) - ((T + x) / x ^ 2) * ( 1 - exp(-x / T)));   
     
    % LF
    tm = abs(pi*H*fmax/P);
    m = -abs(fmax/tm);
    P_LF=- m*k*(T + (-T + tm)*exp(tm/T))*exp(-tm/T);
     
    % EF
    x = abs(2 * P / (2 * H ));
    tm = 2 * exp(1) * fmax / x;
    a = -  exp(1) * fmax / tm;
    b = 1 / tm;
    P_EF = -(a * k ./ (T * b - 1)) * (((T ./ (1 - T * b)) * (exp(-1) - exp(-tm / T))) - tm / exp(1));
   
   result(6,ma)= fmax;
   result(7,ma)= ( P_gen(2)-P_QF)/P_gen(2);
   result(8,ma)= ( P_gen(2)-P_LF)/P_gen(2);
   result(9,ma)= ( P_gen(2)-P_EF)/P_gen(2);
   

   result(7,ma)= P_LF;
   result(8,ma)= P_QF;
   result(9,ma)= P_EF;
   result(10,ma)= P_gen(2);
    
end
result = result';
