function test_suite = test_denoise
initTestSuite;

function test_denoise_default
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h);
  signal_denoised_corr = [0.0741827688375062 0.0791701902526268 0.0760842615272340 0.0750476831774179 0.111279774779568 0.163475053283544 -0.0498263815350539 0.0946073088237311 0.135126562486911 -0.0186090620958193 -0.0748812479991294 -0.103470206059426 0.0234254843251780 0.239772540836257 0.0920583398962312 -0.152180640366891 -0.116682073306156 -0.0459389850762785 -0.00245240039778375 0.0755739164104836 0.102548333512214 0.121099911744184 0.177390507921620 0.240386041553093 0.231105933317157 0.198210924493273 0.175672812990725 0.138822049613034 0.127491615387826 0.121409597186325 0.0994935320130783 0.0760019340865427];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_udwt
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 1);
  signal_denoised_corr = [0.126244615385152 0.0952319712425300 0.0671343607152503 0.0513902979722585 0.0430402732682634 0.0586932575131794 0.0861069751902698 0.0989949047763016 0.0908418658128637 -0.0141454670119059 -0.144791527437026 -0.0185533166035902 0.278351613782131 0.279033706376659 -0.0205012032054263 -0.212367658407976 -0.241484343697995 -0.248582298831059 -0.213374214781743 -0.101963712141109 0.0454248851310567 0.181104333949749 0.275294407293259 0.309076259882059 0.298600450385073 0.259080737796607 0.211123535801718 0.183021783525739 0.171966340866576 0.171616812586097 0.168720006300193 0.151066428184072];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_threshold_low
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 0, [1 3.0 0 0 0 0]);
  signal_denoised_corr = [0.0187742354278351 0.0237616568429558 0.0206757281175629 0.0196391497677469 0.0558712413698966 0.108066519873873 -0.105234914944725 0.0391987754140600 0.0797180290772401 -0.0740175955054904 -0.130289781408801 -0.158878739469097 -0.0319830490844931 0.184364007426586 0.0366498064865601 -0.207589173776562 -0.172090606715827 -0.101347518485950 -0.0578609338074549 0.0201653830008125 0.0471398001025425 0.0656913783345127 0.121981974511949 0.184977508143422 0.175697399907486 0.142802391083602 0.120264279581054 0.0834135162033633 0.0720830819781554 0.0660010637766539 0.0440849986034073 0.0205934006768717];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_thresh_multiplier
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 0, [1 3.5 0 0 0 0]);
  signal_denoised_corr = [0.00563527074803461 0.0110853052404048 0.0101590193471916 0.0116789518546074 0.0354625658443208 0.0691904606426981 -0.0647010252187970 0.0393485097012034 0.0302297746478269 -0.0658230296401878 -0.0947938063374137 -0.147943151851009 -0.0355607514547514 0.143027827800490 0.0126752977970079 -0.200577663821584 -0.149059259007655 -0.0564432101940217 -0.0281365070661950 0.0201021371871464 0.0438412772787373 0.0596866399869512 0.0967101937989458 0.136451641917565 0.130716307107088 0.109146914388131 0.0925200849653435 0.0657607417363412 0.0550584910898860 0.0469636231448182 0.0277268486177313 0.00667135407398081];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_std
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 0, [0 3.0 1 0 0 0]);
  signal_denoised_corr = [0.0686926069658060 0.0706216045196474 0.0719769032529757 0.0743568305131058 0.0754251996534692 0.0763549103855611 0.0783972750744446 0.0807092136475563 0.0763109954998047 0.0693017683604205 0.0628697537191382 0.0547492531677562 0.0755519478401559 0.107931256046656 0.0859959791464885 0.0494376118339224 0.0602059364595448 0.0785077229738383 0.0791999606842265 0.0809410605777517 0.0844652184548917 0.0873749084881920 0.0911535278085727 0.0952027332951270 0.0936316016468421 0.0898878427420561 0.0866734185917041 0.0820709685744921 0.0793481432323076 0.0768306965269240 0.0727995727792393 0.0684196591566048];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_hard
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 0, [0 3.0 0 1 0 0]);
  signal_denoised_corr = [0.0977394160103721 0.0994161560983385 0.0832447407807381 0.0666983311697188 0.177420971595413 0.340230583897110 -0.354597069671295 0.0250017872275015 0.394418485343238 -0.0595745304374512 -0.452401570793399 -0.175707560852101 -0.00622320325130765 0.437867065411816 0.187485346584306 -0.241060664687049 -0.306285896120773 -0.373946536466370 -0.246165924475657 0.00210496326791051 0.0528629966064817 0.0967383656953347 0.275410693617439 0.487298926169970 0.454985253718689 0.348603331393631 0.288205743942248 0.186806596496260 0.172147260405660 0.180050851714681 0.142136445826288 0.104484725401481];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_levels
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 0, [0 3.0 0 0 4 0]);
  signal_denoised_corr = [0.164259992817262 0.156379071218712 0.142212685671703 0.125038963573761 0.150297815252073 0.191536767978636 -0.0381639580765735 0.0881092032192094 0.119629284458486 -0.0406090725365491 -0.105645426731493 -0.141820831994602 -0.0280318977202704 0.173171960129832 0.0117537437282443 -0.247115729957293 -0.206759297285911 -0.123147866042363 -0.0685808245422524 0.0255826360141400 0.0635302930397082 0.0930381970490923 0.165728084463140 0.246884147157615 0.246603211345582 0.220210934934003 0.206436991723089 0.177172675548210 0.178948997433275 0.188010177892750 0.179798128181065 0.170937023676945];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_actual_thresh
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 0, [0 3.0 0 0 0 0.5]);
  signal_denoised_corr = [0.0607099183942295 0.0654351521193524 0.0684154759800610 0.0742018934148454 0.0758845005390013 0.0769511530643110 0.0810856606730252 0.0858023375316036 0.0704706443350518 0.0472060906047587 0.0254329679518446 -0.00154590940405266 0.0598455182579352 0.156556707841878 0.0864272987162393 -0.0287835335280487 0.00606017120154721 0.0659592575432934 0.0713958080495586 0.0812891735076492 0.0953701981347179 0.107554576791239 0.123739146895592 0.141180422640726 0.137085044622601 0.124838366760086 0.114852957437233 0.0997294000571788 0.0922174665178409 0.0857758976557685 0.0737052631031342 0.0605470542090229];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_udwt_threshold_low
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 1, [1 3.0 0 0 0 0]);
  signal_denoised_corr = [0.135039400483741 0.117805175604609 0.0967709584177031 0.0142060292567307 -0.0239840294603812 0.323425861331697 -0.212285200125643 0.166066657685731 0.136653739821785 -0.0361708285655289 -0.244622217319313 -0.0751486112344819 0.279128997196628 0.299915294672821 0.00822389077239383 -0.232180770499244 -0.330137263335199 -0.293955318206172 -0.175538926380835 -0.0733568677543535 0.049241196655251 0.200165899490694 0.304615650610263 0.337325376378116 0.325593984310807 0.282048956150932 0.228861081870546 0.196656880842149 0.180959366486141 0.175210410022406 0.169828050229736 0.155033256209497];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_udwt_thresh_multiplier
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 1, [1 3.5 0 0 0 0]);
  signal_denoised_corr = [0.0479478506866607 0.0160653046305043 -0.012660890293452 -0.0292521383561941 -0.0383355043751224 -0.0239494802109215 0.00200042536526626 0.0135636610003902 0.00399637041195728 -0.100521378500944 -0.229923524965501 -0.102614225576592 0.195850596270724 0.197593413336102 -0.100882406775293 -0.291163630119251 -0.318524834100706 -0.324752887320235 -0.288916218874243 -0.176658530913858 -0.028536592326759 0.108409816572649 0.204063702017061 0.239170248556769 0.230108690684778 0.190119394184444 0.14091827822899 0.11174543739754 0.0991301032767805 0.0977198505254529 0.0937639547688583 0.0745251447941448];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_udwt_std
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 1, [0 3.0 1 0 0 0]);
  signal_denoised_corr = [0.0847626939447046 0.0648669375488877 0.0505127048998841 0.0431477690668965 0.0443458995091662 0.0638361516754724 0.0926698200065443 0.122716357496751 0.135591683864019 0.0377466753027189 -0.0889166586897228 -0.0310700016943258 0.16530654803759 0.237349858169585 0.0577692051497442 -0.137751577705709 -0.18354744395111 -0.188205427540335 -0.157902857480421 -0.055391323576937 0.0791892398460303 0.198068185997372 0.271471422836112 0.282275886815228 0.246689293630916 0.205546705496588 0.16546007731141 0.145130898382968 0.1471329636038 0.142472749823065 0.132163448290946 0.111958195551385];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_udwt_hard
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 1, [0 3.0 0 1 0 0]);
  signal_denoised_corr = [0.135039400483741 0.117805175604609 0.0967709584177031 0.0142060292567307 -0.0239840294603812 0.323425861331697 -0.212285200125643 0.166066657685731 0.136653739821785 -0.0361708285655289 -0.244622217319313 -0.0751486112344819 0.279128997196628 0.299915294672821 0.00822389077239383 -0.232180770499244 -0.330137263335199 -0.293955318206172 -0.175538926380835 -0.0733568677543535 0.049241196655251 0.200165899490694 0.304615650610263 0.337325376378116 0.325593984310807 0.282048956150932 0.228861081870546 0.196656880842149 0.180959366486141 0.175210410022406 0.169828050229736 0.155033256209497];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_udwt_levels
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 1, [0 3.0 0 0 4 0]);
  signal_denoised_corr = [0.137633389000662 0.120676804147327 0.0997827582151432 0.0156985740202669 -0.0251180988153785 0.319788331991522 -0.217919217670089 0.160238201773756 0.131270340429534 -0.0414158027972923 -0.249853610380694 -0.0801267408837784 0.275034335985338 0.296982831400265 0.00620014657281041 -0.234309647934845 -0.33273125185212 -0.296826946748889 -0.178550726178275 -0.0748494125178897 0.0503752660102483 0.203803428830869 0.310249668154709 0.343153832290091 0.330977383703058 0.287293930382695 0.234092474931927 0.201635010491445 0.185054027697432 0.178142873294961 0.171851794429319 0.157162133645098];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);

function test_denoise_udwt_actual_thresh
  signal = makesig('Doppler', 32);
  noise = [1.54421189550395 0.0859311331754255 -1.49159031063761 -0.742301837259857 -1.06158173331999 2.35045722400204 -0.615601881466894 0.748076783703985 -0.192418510588264 0.888610425420721 -0.764849236567874 -1.40226896933876 -1.42237592509150 0.488193909859941 -0.177375156618825 -0.196053487807333 1.41931015064255 0.291584373984183 0.197811053464361 1.58769908997406 -0.804465956349547 0.696624415849607 0.835088165072682 -0.243715140377952 0.215670086403744 -1.16584393148205 -1.14795277889859 0.104874716016494 0.722254032225002 2.58549125261624 -0.666890670701386 0.187331024578940];
  with_noise = signal + noise / 10; 
  h = daubcqf(6);
  [signal_denoised, subtracted_noise, actual_options] = denoise(with_noise, h, 1, [0 3.0 0 0 0 0.5]);
  signal_denoised_corr = [0.126244615385152 0.09523197124253 0.0671343607152503 0.0513902979722585 0.0430402732682634 0.0586932575131794 0.0861069751902698 0.0989949047763016 0.0908418658128637 -0.0141454670119059 -0.144791527437026 -0.0185533166035902 0.278351613782131 0.279033706376659 -0.0205012032054263 -0.212367658407976 -0.241484343697995 -0.248582298831059 -0.213374214781743 -0.101963712141109 0.0454248851310567 0.181104333949749 0.275294407293258 0.309076259882059 0.298600450385073 0.259080737796607 0.211123535801717 0.183021783525739 0.171966340866576 0.171616812586097 0.168720006300193 0.151066428184072];
assertVectorsAlmostEqual(signal_denoised, signal_denoised_corr, 'relative', 0.0001);
