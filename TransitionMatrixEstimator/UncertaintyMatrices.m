%% Solve Systems of Equations and Calculate Upper/Lower/Mean Transition Matrices

%bn = 10; 
SM = zeros(length(PDF2),bn*bn);
for i = 1:length(PDF2) %For every arrival time in BTCs
    IM2 = IM(i,:); %The Probability Contribution Matrix
    b = max(IM2);
    order = (10^floor(log10(b)));
    b = find(IM2 < order);
    IM2(b) = 0;
    PDF2_2 = PDF2(i);
    SM(i,:) = lsqr(IM2, PDF2_2); %Find Least Squares for each row
end

M_guess = zeros(10,10);
M_low= M_guess;
M_high = M_guess;
lim = 1; %constraint, filter out values greater than 1
a = find(SM>lim);
SM1 = SM;
SM1(a) = 0; 
%Top LEFT CORNER
for i = 1:4 %row
   for j = 1:5 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > .075); %filter out data less than .075
       %b = find(a > 0);
       c = find(a < .5); %filter data greater than .5
       b = intersect(b,c); % all elements left after filtering
       dev = std(a(b)); %one stand deviation of remaining elements
       ave = mean(a(b)); %the mean of the remaining elements
       M_guess(i,j) = ave; %we guess the mean
       M_low(i,j) = -dev + ave; % the low matrix is composed of mean minus 1std
       M_high(i,j) = dev + ave; % the high matrix is composed of mean plus 1std
   end    
end

%repeat process for other parts of matrix

%TOP RIGHT CORNER
for i = 1:4 %row
   for j = 7:10 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a <= .125);
       c = find(a > 0);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave;
   end    
end


%Bottom LEFT CORNER
for i = 7:10 %row
   for j = 1:4 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a <= .125);
       c = find(a > 0);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave; 
   end    
end

%Bottom RIGHT CORNER 
for i = 6:10 %row
   for j = 7:10 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > .075);
       %b = find(a > 0);
       c = find(a < .5);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave; 
   end    
end
% 5 and 6

% Elements  5
for i = 1:4 %row
   for j = 5 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > .0);
       c = find(a < .45);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave; 
   end    
end

for i = 5 %row
   for j = 1:5 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > .0);
       c = find(a < .45);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave; 
   end    
end

for i = 5 %row
   for j = 6:10 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > 0);
       c = find(a < .125);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave;
   end    
end

for i = 6:10 %row
   for j = 5 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a < .125);
       c = find(a > 0);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave; 
   end    
end

%Element 6
for i = 1:4 %row
   for j = 6 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > 0);
       c = find(a < .125);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave; 
   end    
end

for i = 6 %row
   for j = 1:5 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > 0);
       c = find(a < .125);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave;
   end    
end

for i = 6 %row
   for j = 6:10 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > .0);
       c = find(a < .45);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = -dev + ave;
       M_high(i,j) = dev + ave; 
   end    
end

for i = 7:10 %row
   for j = 6 %column
       a = SM1(:, (i-1)*10 + j); %gives column of M(i,j) 
       b = find(a > .0);
       c = find(a < .45);
       b = intersect(b,c);
       dev = std(a(b));
       ave = mean(a(b));
       M_guess(i,j) = ave; 
       M_low(i,j) = ave -dev;
       M_high(i,j) = dev + ave;
   end    
end

%% Normalize

%Normalize so the matrix rows sum to 1

% Elements that have high probabilites are given given more weight in
% normalization, meaning they are changed less than elements with low
% probabilities
for i = 1:10
   if sum(M_guess(i,:)) > 1
      weights = 1 - M_guess(i,:); 
      weights = weights/min(weights);
      need = sum(M_guess(i,:))-1;
      for j = 1:10
         M_guess(i,j) = M_guess(i,j) - need*weights(j)/sum(weights); 
      end
   end
end

for i = 1:10
   if sum(M_low(i,:)) > 1
      weights = 1 - M_low(i,:); 
      weights = weights/min(weights);
      need = sum(M_low(i,:))-1;
      for j = 1:10
         M_low(i,j) = M_low(i,j) - need*weights(j)/sum(weights); 
      end
   end
end

for i = 1:10
   if sum(M_high(i,:)) > 1
      weights = 1 - M_high(i,:); 
      weights = weights/min(weights);
      need = sum(M_high(i,:))-1;
      for j = 1:10
         M_high(i,j) = M_high(i,j) - need*weights(j)/sum(weights); 
      end
   end
end

for i = 1:10    
    M_guess(i,:) = M_guess(i,:)/sum(M_guess(i,:));
    M_low(i,:) = M_low(i,:)/sum(M_low(i,:));
    M_high(i,:) = M_high(i,:)/sum(M_high(i,:));
end  

for i = 1:10
    for j = 1:10
        if M_high(i,j) < 0
            M_high(i,j) = 0;
        end
    end
end

           