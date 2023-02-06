//-----------------------------------------------------------------------------
//                              PRÉ-PROCESSAMENTO
//-----------------------------------------------------------------------------
clc
clear
mode(-1)
//-----------------------------------------------------------------------------
//                              FUNÇÕES
//-----------------------------------------------------------------------------
function dxdt = cstr_din(t, x)
    Ca = x(1);
    T = x(2);
    Tj = x(3);
    
    k = alpha*exp(-E/R/T);
    
    dxdt(1) = (F0*Ca0 - F0*Ca - V*k*Ca)/V;                                 // dCa_dt
    dxdt(2) = (rho*Cp*F0*(T0-T) - lambd*V*k*Ca - U*A*(T-Tj))/(rho*Cp*V);   // dT_dt
    dxdt(3) = (rhoj*Cj*Fj*(Tj0-Tj) + U*A*(T-Tj))/(rhoj*Cj*Vj);             // dTj_dt
endfunction

function [f] = cstr_ss(x)
    Ca = x(1);
    T = x(2);
    Tj = x(3);
    
    k = alpha*exp(-E/R/T);
    
    f(1) = F0*(Ca0-Ca)-V*k*Ca;                              // Ca
    f(2) = rho*Cp*F0*(T0-T)-lambd*V*k*Ca-U*A*(T-Tj);        // T
    f(3) = rhoj*Cj*Fj*(Tj0-Tj)+U*A*(T-Tj);                  // Tj
endfunction
//-----------------------------------------------------------------------------
//                              INICIAÇÃO DE VALORES
//-----------------------------------------------------------------------------
// Valores Parâmetros do CSTR
F0 = 40;                 // ft^3/hr
F = 40;                  // ft^3/hr
V = 48;                  // ft^3
Fj = 49.9;               // ft^3/hr
Vj = 3.85;               // ft^3
alpha = 7.08e10;         // hr^-1
E = 30000;               // BTU/mol
R = 1.99;                // BTU/mol ºR
U = 150;                 // BTU/hr-ft^2-ºR
A = 250;                 // ft^2
lambd = -30000;          // BTU/mol
Cp = 0.75;               // BTU/lbm-ºR
Cj = 1;                  // BTU/lbm-ºR
rho = 50;                // lmb/ft^3
rhoj = 62.3;             // lmb/ft^3
// Condições iniciais
Ca0 = 0.50;              // mol/ft^3
Tj0 = 530;               // ºR
T0 = 530;                // ºR
//-----------------------------------------------------------------------------
//                              PROGRAMA PRINCIPAL
//-----------------------------------------------------------------------------
       
t = [0:0.1:10];                 // vetor de tempo
CI = [Ca0; T0; Tj0];            // vetor de condição inicial

//ESTADO ESTACIONÁRIO
[sol, fsol, info] = fsolve(CI, cstr_ss)
if info==1 then
    disp('Solução:')
    disp(sol)
else
    disp('tente novamente')
end

//DINÂMICO
y = ode(CI, t(1), t, cstr_din)
disp(y)

subplot(3,1,1),plot(t,y(1,:))
subplot(3,1,2),plot(t,y(2,:))
subplot(3,1,3),plot(t,y(3,:))
