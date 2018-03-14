% sekcja czyszczenia
clear all; clc; close all;

% nazwy obrazów w sekwencji 
names = {'IM17.jpg', 'IM18.jpg', 'IM19.jpg', 'IM20.jpg', 'IM21.jpg', 'IM22.jpg'};

% rozmiary komór, ptrzebne do funkcji bwareaopen
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
 % trójwymiarowa macierz zawieraj¹ca ca³¹ sekwencjê
 ALL = [];
 
 % macierz rozmiarów rzeczywistego obrazu i wysegmentowanych komór
 pixels = [];

 for i = 1 : 1 : length(names)
 
    % wczytanie obrazu
    image = dicomread(names{i});
    I(:, :) = double(image);

    % obliczenie rozmiarów obrazu
    imagesize = length(image(1, :)) * length(image(:, 1)); % mierzymy 1. kolumnê i 1. wiersz
    
    % obliczenie pixeli t³a
    zeropixels = size(find(I == 0));
    
    % rzeczywisty rozmiat obrazu (bez t³a)
    imagepixels = imagesize - zeropixels(1);
    
    % wyœwietlenie oryginalnego obrazu
    figure();
    imshow(I, []); 
    title('Original image');
    hold on;

    data = reshape(I, [], 1);

    % zapisanie wspó³rzêdnych komór (punkt g³ówny)- wy³¹cznie dla
    % pierwszego obrazu
    [xo, yo] = ginput(1);
    xo = round(xo);
    yo = round(yo);

    % mo¿liwe wartoœci komór po lustrzanym odbiciu
    xo2 = 512 - xo;
    % yo pozostaje bez zmian

    % naniesienie punktu komór
    hold on; 
    p = plot(xo, yo, '.');
    set(p, 'MarkerSize', 15);
    set(p, 'Color', 'g');

    % naniesienie punktu komór po symetrycznym odbiciu
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
    pause(0.2); % potrzebne do wyœwitetlenia kropki
    hold on;

    % nie zmieniaæ
    klasa_obszaru = 1;
    klasa_pozaobszarem = 1;
    ilosc_klas = 2;

    % proces segmentacji dla wskazanej czêœci obrazu
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

    % erozja (wyodrêbienie komór jako jednego obszaru)
    erode1 = imerode(BW1, se);

%     % wyœwietlenie
%     figure();
%     imshow(erode1);
%     title('Isolated brain ventricles');
%     hold on;

    % wyciêcie wszystkiego oprócz komór
    
    vent1 = bwareaopen(erode1, sizes(i));  % drugi argument mówi o tym jak du¿y element 
                                  % powinien byæ zostawiamy

    % powtórzenie procesu dla czêœci po drugiej stronie osi symetrii
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

    % zsumowanie lustrzanych odbiæ
    BW = BW1 + BW2;

    % erozja (wyodrêbienie komór jako jednego obszaru)
    erode = imerode(BW, se);

%     % wyœwietlenie
%     figure();
%     imshow(erode);
%     title('Isolated brain ventricles');
%     hold on;

    % wyciêcie wszystkiego oprócz komór
    vent = bwareaopen(erode, sizes(i)); % drugi argument mówi o tym jak du¿y element 
                                        % powinien byæ zostawiamy

    % odtworzenie oryginalnych rozmiarów komór
    dilate = imdilate(vent, se);

    % wyliczenie liczby pixeli
    blackpixels = sum(sum(dilate == 0));
    whitepixels = sum(sum(dilate));
    
    whiteproc = (whitepixels/(whitepixels+blackpixels))*100; % procentowo
    blackproc = (blackpixels/(whitepixels+blackpixels))*100; % procentowo

    % wyœwietlenie wy³¹cznie komór
    figure();
    imshow(dilate);
    title('Only brain ventricles');
    xlabel(['White pixels: ', num2str(whiteproc), ' %']);
    ylabel(['Black pixels: ', num2str(blackproc), ' %']);
    
    % zapisanie wartoœci lokalnych do globalnych zmiennych
    ALL(:, :, i) = dilate;
    pixels(i, 1) = imagepixels();
    pixels(i, 2) = whitepixels;
    
end % pêtla

% obliczenie objêtoœci g³owy
head = 0;
for k = 1 : length(pixels(:, 1))
    head = head + pixels(k, 1);
end

%obliczenie objêtoœci wystegmentowanych komór
vent = 0;
for k = 1 : length(pixels(:, 2))
    vent = vent + pixels(k, 2);
end

% obliczenie objêtoœci wysegmentowanych komór
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
