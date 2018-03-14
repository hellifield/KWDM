% sekcja czyszczenia
clear all; clc; close all;

% nazwy obraz�w w sekwencji 
names = {'IM17.jpg', 'IM18.jpg', 'IM19.jpg', 'IM20.jpg', 'IM21.jpg', 'IM22.jpg'};

% rozmiary kom�r, ptrzebne do funkcji bwareaopen
sizes = [20 100 1500 1500 500 500];
%sizes = [240 100 1500 1500 500 500 100 100 100 100];

% element strukturalny
se = [1 1 1 1 1 1 1 1;
      1 0 0 0 0 0 0 1;
      1 0 0 0 0 0 0 1;
      1 0 0 0 0 0 0 1;
      1 0 0 0 0 0 0 1;
      1 0 0 0 0 0 0 1;
      1 1 1 1 1 1 1 1];
 
 % Inicjalizacja 
 % tr�jwymiarowa macierz zawieraj�ca ca�� sekwencj�
 ALL = [];
 
 % macierz rozmiar�w rzeczywistego obrazu i wysegmentowanych kom�r
 pixels = [];

 for i = 1 : 1 : length(names)
 
    % wczytanie obrazu
    image = dicomread(names{i});
    I(:, :) = double(image);

    % obliczenie rozmiar�w obrazu
    imagesize = length(image(1, :)) * length(image(:, 1)); % mierzymy 1. kolumn� i 1. wiersz
    
    % obliczenie pixeli t�a
    zeropixels = size(find(I == 0));
    
    % rzeczywisty rozmiat obrazu (bez t�a)
    imagepixels = imagesize - zeropixels(1);
    
    % wy�wietlenie oryginalnego obrazu
    figure();
    imshow(I, []); 
    title('Original image');
    hold on;

    data = reshape(I, [], 1);

    % zapisanie wsp�rz�dnych kom�r (punkt g��wny)- wy��cznie dla
    % pierwszego obrazu
    [xo, yo] = ginput(1);
    xo = round(xo);
    yo = round(yo);

    % mo�liwe warto�ci kom�r po lustrzanym odbiciu
    xo2 = 512 - xo;
    % yo pozostaje bez zmian

    % naniesienie punktu kom�r
    hold on; 
    p = plot(xo, yo, '.');
    set(p, 'MarkerSize', 15);
    set(p, 'Color', 'g');

    % naniesienie punktu kom�r po symetrycznym odbiciu
    hold on; 
    p = plot(xo2, yo, '.');
    set(p, 'MarkerSize', 15);
    set(p, 'Color', 'g');

    % wskazanie punktu odniesienia
    [xp, yp] = ginput(1);
    xp = round(xp);
    yp = round(yp);

    % naniesienie punktu odniesienia
    hold on; 
    p = plot(xp, yp, '.');
    set(p, 'MarkerSize', 15);
    set(p, 'Color', 'r');
    pause(0.2); % potrzebne do wy�witetlenia kropki
    hold on;

    % nie zmienia�
    klasa_obszaru = 1;
    klasa_pozaobszarem = 1;
    ilosc_klas = 2;

    % proces segmentacji dla wskazanej cz�ci obrazu
    while (klasa_obszaru == klasa_pozaobszarem),
        [center, U, obj_fcn] = fcm(data, ilosc_klas, [NaN NaN NaN 0]);
        [maxU, nrklasy] = max(U', [], 2);
        I_klasy = zeros(size(I));
        I_nr_klas = zeros(size(I));
        maxnrklasy = max(nrklasy);

        for j = 1 : maxnrklasy,
            idx = find(nrklasy == j);
            I_klasy(idx) = (ones(length(idx), 1)*center(j));
            I_nr_klas(idx) = j;
        end

        klasa_obszaru = I_nr_klas(yo, xo);
        klasa_pozaobszarem = I_nr_klas(yp, xp);
        ilosc_klas = ilosc_klas + 1;

        if ilosc_klas > 5 
            break;
        end
    end

    idx = find(I_nr_klas == klasa_obszaru);

    macierz = zeros(size(I));
    macierz(idx) = 1;

    xp1 = xo;
    yp1 = yo;

    xp2 = xp1;
    yp2 = yp1;

    wiersz = [xp1 xp2];
    kolumna = [yp1 yp2];
    czesc_klasy1 = bwselect(macierz, wiersz, kolumna, 8);

    BW1 = bwmorph(czesc_klasy1, 'clean', 5);

    % erozja (wyodr�bienie kom�r jako jednego obszaru)
    erode1 = imerode(BW1, se);

%     % wy�wietlenie
%     figure();
%     imshow(erode1);
%     title('Isolated brain ventricles');
%     hold on;

    % wyci�cie wszystkiego opr�cz kom�r
    
    vent1 = bwareaopen(erode1, sizes(i));  % drugi argument m�wi o tym jak du�y element 
                                  % powinien by� zostawiamy

    % powt�rzenie procesu dla cz�ci po drugiej stronie osi symetrii
    xo = xo2; % zmiana punktu wskazanej klasy

    while (klasa_obszaru == klasa_pozaobszarem),
        [center, U, obj_fcn] = fcm(data, ilosc_klas, [NaN NaN NaN 0]);
        [maxU, nrklasy] = max(U', [], 2);
        I_klasy = zeros(size(I));
        I_nr_klas = zeros(size(I));
        maxnrklasy = max(nrklasy);

        for j = 1 : maxnrklasy,
            idx = find(nrklasy == j);
            I_klasy(idx) = (ones(length(idx), 1)*center(j));
            I_nr_klas(idx) = j;
        end

        klasa_obszaru = I_nr_klas(yo, xo);
        klasa_pozaobszarem = I_nr_klas(yp, xp);
        ilosc_klas = ilosc_klas + 1;

        if ilosc_klas > 5 
            break;
        end
    end

    idx = find(I_nr_klas == klasa_obszaru);

    macierz = zeros(size(I));
    macierz(idx) = 1;

    xp1 = xo;
    yp1 = yo;

    xp2 = xp1;
    yp2 = yp1;

    wiersz = [xp1 xp2];
    kolumna = [yp1 yp2];
    czesc_klasy2 = bwselect(macierz, wiersz, kolumna, 8);

    BW2 = bwmorph(czesc_klasy2, 'clean', 5);

    % zsumowanie lustrzanych odbi�
    BW = BW1 + BW2;

    % erozja (wyodr�bienie kom�r jako jednego obszaru)
    erode = imerode(BW, se);

%     % wy�wietlenie
%     figure();
%     imshow(erode);
%     title('Isolated brain ventricles');
%     hold on;

    % wyci�cie wszystkiego opr�cz kom�r
    vent = bwareaopen(erode, sizes(i)); % drugi argument m�wi o tym jak du�y element 
                                        % powinien by� zostawiamy

    % odtworzenie oryginalnych rozmiar�w kom�r
    dilate = imdilate(vent, se);

    % wyliczenie liczby pixeli
    blackpixels = sum(sum(dilate == 0));
    whitepixels = sum(sum(dilate));
    
    whiteproc = (whitepixels/(whitepixels+blackpixels))*100; % procentowo
    blackproc = (blackpixels/(whitepixels+blackpixels))*100; % procentowo

    % wy�wietlenie wy��cznie kom�r
    figure();
    imshow(dilate);
    title('Only brain ventricles');
    xlabel(['White pixels: ', num2str(whiteproc), ' %']);
    ylabel(['Black pixels: ', num2str(blackproc), ' %']);
    
    % zapisanie warto�ci lokalnych do globalnych zmiennych
    ALL(:, :, i) = dilate;
    pixels(i, 1) = imagepixels();
    pixels(i, 2) = whitepixels;
    
end % p�tla

% obliczenie obj�to�ci g�owy
head = 0;
for k = 1 : length(pixels(:, 1))
    head = head + pixels(k, 1);
end

%obliczenie obj�to�ci wystegmentowanych kom�r
vent = 0;
for k = 1 : length(pixels(:, 2))
    vent = vent + pixels(k, 2);
end

% obliczenie obj�to�ci wysegmentowanych kom�r
volume = (vent/head)*100; % procentowo

%% rekonstrukcja 3D
figure();
D = squeeze(ALL);
h = vol3d('cdata',D,'texture','3D');
view(3);  
axis tight;  
daspect([1 1 0.4])
alphamap('rampup');
alphamap(0.06 .* alphamap);
ylabel(['Ventricles volume: ', num2str(volume), ' %']);
