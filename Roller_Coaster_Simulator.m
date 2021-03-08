%Roller Coaster Simualtor

%Authors: 

%Erick Alberto Bustos Cruz A01378966
%Eduardo Rodríguez López A01749381
%Carlos Antonio Pazos Reyes A01378262
% Nancy Lesly García Jiménez A01378043
% Jared Abraham Flores Guarneros A01379868

function Proyecto
    fig = uifigure("Name","Montaña Rusa");   %Crea el display
    ax = uiaxes("Parent",fig,"Units","pixels","Position",[5,5,400,350],"ylim",[0,45],"xlim",[0,55]);   %Crea el eje de grafica en el display
    ax_2 = uiaxes("Parent",fig,"Units","pixels","Position",[410,5,130,350],"ylim",[0,10000],"xlim",[0,10]); % Crea eje barras en display
    vel_label = uilabel(fig,"Position",[5,360,70,20],"Text","Vel. Inicial: ");  % Crea el texto de la velocidad en display
        vel = uieditfield(fig,"numeric","Position",[65,360,60,20]);    % Crea cuadro editable vel display
        vel.Value = 20;      %Pone un valor por default para velocidad
    m_mass_label = uilabel(fig,"Position",[445,360,60,20],"Text","Masa: ");    % Crea texto de masa en display
        m_mass = uieditfield(fig,"numeric","Position",[480,360,60,20]);       % Crea cuadro editable masa en display
        m_mass.Value = 20;           % Pone valor default masa
    coef_fr_label = uilabel(fig,"Position",[240,360,120,20],"Text","Coeficiente Friccion: ");  %Crea texto de coeficiente friccion en display
        coef_fr = uieditfield(fig,"numeric","Position",[350,360,60,20]);      %Crea cuadro editable para coeficiente friccion en display
        coef_fr.Value = 0.1;     %Pone valor default para coeficiente friccion
    coef_aire_label = uilabel(fig,"Position",[145,360,60,20],"Text","k: ");  %Crea texto de coeficiente aire en display
        coef_aire = uieditfield(fig,"numeric","Position",[155,360,60,20]);  %Crea cuadro editable para coeficiente aire en display
        coef_aire.Value = 0.1;     %Pone valor default para coeficiente aire
    btn = uibutton(fig,"Push","Text","Plot","Position",[225,390,100,20],... % Crea boton con texto plot en display
            "ButtonPushedFcn",@(btn,event)ploteoButtonPushed(ax,ax_2,coef_fr,coef_aire,vel,m_mass));  %Alista la funcion que realiza al push button
    
end

function ploteoButtonPushed(ax,ax_2,coef_fr,coef_aire,vel,m_mass)
    %--------------------------- FUNCIONES------------------------------------
    %PRIMER TROZO
    syms f(x)  % Se crea una funcion simbólica de la pista 
    f(x) = -(127 ./ 625) .* x .^2 + 2.8509 .* x +30;
    dfx = diff(f,x);% Se deriva la función para su futura implementación
    d2fx = diff(dfx,x);% Se encuentra la segunda derivada

    %SEGUNDO TROZO
    syms h(x)
    h(x) = (9380.53 .* sin(0.01 .* (x))) - (1555.94 .* cos(0.01 .* (x))) - ...
            (4603.35 .* sin(2 .* 0.01 .* (x))) + (1596.73 .* cos(2 .* 0.01 .* (x)));
    dhx = diff(h,x); 
    d2hx = diff(dhx,x); 

    %TERCER TROZO
    syms g(x)
    g(x) = -0.0021 .* x .^ 3 + 0.151 .* x .^2 - 0.3888 .* x - 70.695;
    dgx = diff(g,x);
    d2gx = diff(dgx,x);

    %-------------------------VARIABLES INICIALES------------------------------

    gravedad = 9.8;               % Gravedad
    delta_t = 0.1;         % Tamaño de paso
    v = vel.Value;                 % Velocidad inicial
    xprima = 0;            % x en eje primo
    xpos = 0;              % x en eje normal
    ypos = double(f(xpos));% y en el eje normal
    masa = m_mass.Value;              % Masa del objeto
    muK = coef_fr.Value;             % Coeficiente de friccion cinetica
    coef_resistencia = coef_aire.Value; % Coeficiente de resistencia del aire
    r = ((1 + (double(dfx(xpos))) .^ 2).^(3./2)) ./ (double(d2fx(xpos))); %Radio
    angulo = atan(double(dfx(xpos))); %angulo
    peso = -masa .* gravedad .* sin(angulo); % Fuerza referente al peso
    friccion = -muK .* masa .*( gravedad .* cos(angulo) +  ((v.^2)./r));
    resistencia_aire = - coef_resistencia .* (v .^2);
    sumadefuerzasenxprima = peso + friccion + resistencia_aire; % Suma de fuerzas inicial 
    aprima = (sumadefuerzasenxprima) ./ masa; % Aceleración inicial


    %--------------------------PLOTEAR PISTA-----------------------------------
    % f(x)
    parte1x = 0:0.1:12; % Vector de valores en x de la primera parte de la pista
    parte1y = ((-127 ./ 625).* (parte1x .^ 2)) + (2.8509 .* parte1x) + 30; % Vector de valores en y

    % h(x)
    parte2x = 12:0.1:35;
    parte2y = (9380.53 .* sin(0.01 .* parte2x)) - (1555.94 .* cos(0.01 .* parte2x)) - ...
            (4603.35 .* sin(2 .* 0.01 .* parte2x)) + (1596.73 .* cos(2 .* 0.01 .* parte2x));
    % g(x)
    parte3x = 35:0.1:54;
    parte3y = -0.0021 .* (parte3x .^3) + 0.151 .* (parte3x .^2) - 0.3888 .* parte3x - 70.695 ; 


    %------------------------INICIO DE LA SIMULACION---------------------------

    while xpos <= 54 && xpos >= 0
        % Simulación para f(x)
        if (xpos >= 0) && (12 >= xpos)
       
            % Se establecen valores "-1"
            xprimaanterior = xprima;
            xposanterior = xpos;
            yposanterior = ypos;
            vanterior = v;
            aprimaanterior = aprima;
            ranterior = r;

            % Se calcula xprima siguiente usando Velocity Verlet
            xprima = xprimaanterior + vanterior .* delta_t + 0.5 .* aprimaanterior .* (delta_t .^2);

            % Se calculan los nuevos valores de xpos y ypos
            xpos = xposanterior + (xprima - xprimaanterior) .* cos(angulo);
            ypos = double(f(xpos));
            
            % Se dibujan los puntos haciendo uso de la animatedline
            plot(ax,parte1x,parte1y,"k",parte2x,parte2y,"k",parte3x, parte3y,"k",xpos,ypos,"o")
            E_pot = masa*gravedad*ypos;
            E_cin = 0.5*masa*(vanterior.^2);
            E_total = E_pot+E_cin;
            fuerzas = [E_pot;E_cin;E_total];
            c = categorical(["Potencial","Cinetica","Mecanica"]);
            b = bar(ax_2,c,fuerzas,"FaceColor","flat");
            b.CData(2,:) = [0 1 0];
            b.CData(3,:) = [1 0 0];
            pause(0)
            
            % Se calculan nuevos valores de radio, velocidad y angulo
            r = ((1 + (double(dfx(xpos))) .^ 2).^(3./2)) ./ (double(d2fx(xpos)));
            v = sqrt(((vanterior .^2) .*(1- sign(vanterior) .* ((muK ./ ranterior) + (coef_resistencia ./ masa)).* ...
                (xprima - xprimaanterior)) - 2 .* gravedad .* (ypos - yposanterior) - 2.* sign(vanterior) .* gravedad .* muK .* ...
                (xpos -xposanterior))./ (1 + sign(vanterior) .*((muK ./ r) + (coef_resistencia ./ masa)) .* (xprima - xprimaanterior)));
            angulo = atan(double(dfx(xpos))); 
            
            % Se calcula el valor adentro de raiz para poder hacer que el
            % carrito cambie de dirección
            adentro_raiz = ((vanterior .^2) .*(1-((muK ./ ranterior) + (coef_resistencia ./ masa)).* ...
                (xprima - xprimaanterior)) - 2 .* gravedad .* (ypos - yposanterior) - 2 .* gravedad .* muK .* ...
                (xpos -xposanterior))./ (1 + ((muK ./ r) + (coef_resistencia ./ masa)) .* (xprima - xprimaanterior));

            if adentro_raiz <= 0.2
                v = -vanterior;
                xprima = xprimaanterior;
                xpos = xposanterior;
                ypos = yposanterior;
                aprima = aprimaanterior; 
            else
                v = sign(vanterior) .* v;
                if sign(vanterior) == 1
                    peso = masa .* gravedad .* sin(angulo); % Fuerza referente al peso
                    friccion = muK .* masa .*( gravedad .* cos(angulo) +  ((v.^2)./r));
                    resistencia_aire = coef_resistencia .* (v .^2);
                    sumadefuerzasenxprima = - peso - friccion - resistencia_aire; % Suma de fuerzas inicial 
                    aprima = (sumadefuerzasenxprima) ./ masa; % Aceleración inicial
                else
                    peso = masa .* gravedad .* sin(angulo); % Fuerza referente al peso
                    friccion = muK .* masa .*( gravedad .* cos(angulo) +  ((v.^2)./r));
                    resistencia_aire =  coef_resistencia .* (v .^2);
                    sumadefuerzasenxprima = peso + friccion + resistencia_aire; % Suma de fuerzas inicial 
                    aprima = (sumadefuerzasenxprima) ./ masa; % Aceleración inicial
                end
            end
        
    end
    %-------------------------------------------------------------------------------------------
        if (xpos > 12) && (35 >= xpos)
            % Se establecen valores "-1"
            xprimaanterior = xprima;
            xposanterior = xpos;
            yposanterior = ypos;
            vanterior = v;
            aprimaanterior = aprima;
            ranterior = r;

            % Se calcula xprima siguiente usando Velocity Verlet
            xprima = xprimaanterior + vanterior .* delta_t + 0.5 .* aprimaanterior .* (delta_t .^2);

            % Se calculan los nuevos valores de xpos y ypos
            xpos = xposanterior + (xprima - xprimaanterior) .* cos(angulo);
            ypos = double(h(xpos));

            % Se dibujan los puntos haciendo uso de la animatedline
            plot(ax,parte1x,parte1y,"k",parte2x,parte2y,"k",parte3x, parte3y,"k",xpos,ypos,"o")
            E_pot = masa*gravedad*ypos;
            E_cin = 0.5*masa*(vanterior.^2);
            E_total = E_pot+E_cin;
            fuerzas = [E_pot;E_cin;E_total];
            c = categorical(["Potencial","Cinetica","Mecanica"]);
            b = bar(ax_2,c,fuerzas,"FaceColor","flat");
            b.CData(2,:) = [0 1 0];
            b.CData(3,:) = [1 0 0];
            pause(0)
        
            % Se calculan nuevos valores de radio, velocidad y angulo
            r = ((1 + (double(dhx(xpos))) .^ 2).^(3./2)) ./ (double(d2hx(xpos)));
            v = sqrt(((vanterior .^2) .*(1- sign(vanterior) .* ((muK ./ ranterior) + (coef_resistencia ./ masa)).* ...
                (xprima - xprimaanterior)) - 2 .* gravedad .* (ypos - yposanterior) - 2.* sign(vanterior) .* gravedad .* muK .* ...
                (xpos -xposanterior))./ (1 + sign(vanterior) .*((muK ./ r) + (coef_resistencia ./ masa)) .* (xprima - xprimaanterior)));
            angulo = atan(double(dhx(xpos))); 
            % Se calcula el valor adentro de raiz para poder hacer que el
            % carrito cambie de dirección
            adentro_raiz = ((vanterior .^2) .*(1-((muK ./ ranterior) + (coef_resistencia ./ masa)).* ...
                (xprima - xprimaanterior)) - 2 .* gravedad .* (ypos - yposanterior) - 2 .* gravedad .* muK .* ...
                (xpos -xposanterior))./ (1 + ((muK ./ r) + (coef_resistencia ./ masa)) .* (xprima - xprimaanterior));

            if adentro_raiz <= 0.2
                v = -vanterior;
                xprima = xprimaanterior;
                xpos = xposanterior;
                ypos = yposanterior;
                aprima = aprimaanterior; 
            else
                v = sign(vanterior) .* v;
                if sign(vanterior) == 1
                    peso = masa .* gravedad .* sin(angulo); % Fuerza referente al peso
                    friccion = muK .* masa .*( gravedad .* cos(angulo) +  ((v.^2)./r));
                    resistencia_aire = coef_resistencia .* (v .^2);
                    sumadefuerzasenxprima = -peso - friccion - resistencia_aire; % Suma de fuerzas inicial 
                    aprima = (sumadefuerzasenxprima) ./ masa; % Aceleración inicial
                else
                    peso = - masa .* gravedad .* sin(angulo); % Fuerza referente al peso
                    friccion = muK .* masa .*( gravedad .* cos(angulo) -  ((v.^2)./r));
                    resistencia_aire =  coef_resistencia .* (v .^2);
                    sumadefuerzasenxprima = peso + friccion + resistencia_aire; % Suma de fuerzas inicial 
                    aprima = (sumadefuerzasenxprima) ./ masa; % Aceleración inicial
                end
            end
        end
   %-----------------------------------------------------------------------
        if (xpos > 35) && (54 >= xpos)
            % Se establecen valores "-1"
            xprimaanterior = xprima;
            xposanterior = xpos;
            yposanterior = ypos;
            vanterior = v;
            aprimaanterior = aprima;
            ranterior = r;

            % Se calcula xprima siguiente usando Velocity Verlet
            xprima = xprimaanterior + vanterior .* delta_t + 0.5 .* aprimaanterior .* (delta_t .^2);

            % Se calculan los nuevos valores de xpos y ypos
            xpos = xposanterior + (xprima - xprimaanterior) .* cos(angulo);
            ypos = double(g(xpos));

            % Se dibujan los puntos haciendo uso de la animatedline
            plot(ax,parte1x,parte1y,"k",parte2x,parte2y,"k",parte3x, parte3y,"k",xpos,ypos,"o")
            E_pot = masa*gravedad*ypos;
            E_cin = 0.5*masa*(vanterior.^2);
            E_total = E_pot+E_cin;
            fuerzas = [E_pot;E_cin;E_total];
            c = categorical(["Potencial","Cinetica","Mecanica"]);
            b = bar(ax_2,c,fuerzas,"FaceColor","flat");
            b.CData(2,:) = [0 1 0];
            b.CData(3,:) = [1 0 0];
            pause(0)
        
            % Se calculan nuevos valores de radio, velocidad y angulo
            r = ((1 + (double(dgx(xpos))) .^ 2).^(3./2)) ./ (double(d2gx(xpos)));
            v = sqrt(((vanterior .^2) .*(1- sign(vanterior) .* ((muK ./ ranterior) + (coef_resistencia ./ masa)).* ...
                (xprima - xprimaanterior)) - 2 .* gravedad .* (ypos - yposanterior) - 2.* sign(vanterior) .* gravedad .* muK .* ...
                (xpos -xposanterior))./ (1 + sign(vanterior) .*((muK ./ r) + (coef_resistencia ./ masa)) .* (xprima - xprimaanterior)));
            angulo = atan(double(dgx(xpos))); 
            % Se calcula el valor adentro de raiz para poder hacer que el
            % carrito cambie de dirección
            adentro_raiz = ((vanterior .^2) .*(1-((muK ./ ranterior) + (coef_resistencia ./ masa)).* ...
                (xprima - xprimaanterior)) - 2 .* gravedad .* (ypos - yposanterior) - 2 .* gravedad .* muK .* ...
                (xpos -xposanterior))./ (1 + ((muK ./ r) + (coef_resistencia ./ masa)) .* (xprima - xprimaanterior));

            if adentro_raiz <= 0.2
                v = -vanterior;
                xprima = xprimaanterior;
                xpos = xposanterior;
                ypos = yposanterior;
                aprima = aprimaanterior; 
            else
                v = sign(vanterior) .* v;
                if sign(vanterior) == 1
                    peso = masa .* gravedad .* sin(angulo); % Fuerza referente al peso
                    friccion = muK .* masa .*( gravedad .* cos(angulo) +  ((v.^2)./r));
                    resistencia_aire = coef_resistencia .* (v .^2);
                    sumadefuerzasenxprima = -peso - friccion - resistencia_aire; % Suma de fuerzas inicial 
                    aprima = (sumadefuerzasenxprima) ./ masa; % Aceleración inicial
                else
                    peso = -masa .* gravedad .* sin(angulo); % Fuerza referente al peso
                    friccion = muK .* masa .*( gravedad .* cos(angulo) -  ((v.^2)./r));
                    resistencia_aire =  coef_resistencia .* (v .^2);
                    sumadefuerzasenxprima = -peso - friccion - resistencia_aire; % Suma de fuerzas inicial 
                    aprima = (sumadefuerzasenxprima) ./ masa; % Aceleración inicial
                end
            end
        end
    end
end




