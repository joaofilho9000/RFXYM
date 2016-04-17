!  Copyright 2014 João Batista dos Santos Filho <joao@jbsantosfilho.com>
! programa campo aleatorio na direção z
program RFXYM
    
    implicit none
    save                                                                                                                         
    
    type :: Spin
        double precision, dimension(1:3) :: s
    end type Spin

    type :: CampoAleatorio
        double precision, dimension(1:3) :: h
    end type CampoAleatorio
    
    type :: Quantidades
        double precision,dimension(1:3) :: magnetizacao !magnetização instantanea
        double precision :: energia !Energia instantanea
        double precision :: helicidade1
        double precision :: helicidade2
        double precision :: helicidade3
    end type Quantidades

    type :: Configuracao
        integer :: l
        integer :: numeroCoordenacao
        integer :: numeroSitios
        integer :: numeroTermalizacao
        integer :: numeroPassosMC
        integer :: numeroAmostras
        integer :: numeroIntervalosTeperatura
        double precision :: tInicial
        double precision :: tFinal
        double precision :: tsimulacao
        double precision :: deltaT
        double precision :: concentracao 
        double precision :: alfa
        integer :: histograma
        double precision :: hInicial
        double precision :: hFinal
    end type Configuracao

    !decalração de variaveis 
    intrinsic random_seed, random_number
    
    type(Spin), dimension(:), allocatable :: rede
    type(CampoAleatorio), dimension(:), allocatable :: hr
    type(Quantidades) :: quantidadesTermodinamicas
    type(Configuracao) :: sistema 
    integer, dimension(:,:), allocatable :: vizinhos
    integer, dimension(:,:), allocatable :: segvizinhos
    integer, dimension(:,:), allocatable :: ligacao
    real :: rando



    integer :: i, n
    real :: t
    real :: tempoInicial, tempoFinal, tempoCPU
    real :: tempoInicialTotal, tempoFinalTotal, tempoCPUTotal

  !---------------------------------------------------------------------
    integer :: magnetizacao, eneJ1, eneJ2
    double precision :: energia
    double precision :: somaMag, somaMag2, somaMag4
    double precision :: somaEne, somaEne2, somaEne4
    double precision :: mediaMag,  mediaMag2, mediaMag4
    double precision :: mediaEne,mediaEne2, mediaEne4
    double precision :: susceptibilidade, cumuM
    double precision :: calor
    double precision::t0, tin,tfi,dt
  !---------------------------------------------------------------------
    call cpu_time(tempoInicialTotal )
    open(1,file = 'dados.dat')
  !  call lerDados(sistema)
    call leargumentos(sistema)

    !Abrindo arquivos
    open(2,file = 'hist.dat')
    open(3,file = 'tempo.dat')
   
    !inicilização das variaveis 
    allocate(rede(1:sistema%numeroSitios))
    allocate(hr(1:sistema%numeroSitios))
    allocate(vizinhos(1:sistema%numeroCoordenacao,1:sistema%numeroSitios))
    allocate(ligacao(1:sistema%numeroCoordenacao/2,1:sistema%numeroSitios))
    forall(i = 1:sistema%numeroSitios)
        Rede(i)%S=(/1,0,0/)
    end forall
   
    call init_random_seed()   
    ! call random_seed(PUT=seed)
     
    call marcarVizinhos
    call cpu_time(tempoInicial)
    !repecião nas amostras
    T = sistema%tSimulacao
    if (sistema%histograma==1) then
        do n=1, sistema%numeroAmostras
            write(*,*) "simulando amostra " , n
            call IniciaCampoAleatorio !gera amostra
            call inicializar ! calcula energia e magnetização inicial
            !      repetições para equilibrio
            do i = 1, sistema%numeroTermalizacao
                call metropolis(T)
                !call wolff
                call superrelaxacao
            end do
            !repetições para media
            do i = 1, sistema%numeroPassosMC
                call metropolis(T)
                !call wolff
                call superrelaxacao
                call gravarHistograma
            end do

            call cpu_time(tempoFinal)
            TempoCPU=tempoFinal-tempoInicial
            write(3,*) n, tempoCPU
        end do 
        !verificando o tempo gasto
        call cpu_time(tempoFinalTotal)
        TempoCPUTotal = tempoFinalTotal-tempoInicialTotal
        !fechando arquivos
        do i = 1, 3
            close(i)
        end do
     else
        t0 = sistema%tInicial
        tfi = sistema%Tfinal
        dt= sistema%deltaT
        write(*,*) '#iniciando simulação com :'
        write(*,*) '# L =', sistema%l
        write(*,*) '# p =', sistema%alfa
        write(*,*) '# Tin =', Tin
        write(*,*) '# Tfi =', Tfi
        write(*,*) '# dt =', dt
        write(*,*) '# h =', sistema%alfa
!        write(*,*) '# J2 =', J2
        write(*,*) '# numero Termalizacao =', sistema%numeroTermalizacao
        write(*,*) '# numero Passos MC =', sistema%numeroPassosMC
!        write(*,*) '# histograma =', histograma
        call IniciaCampoAleatorio !gera amostra
        call inicializar

        do while (t0<=tfi)
            call iniciaVariaveis
            !repetiÃ§Ãµes termalizaÃ§Ã£o
            do i = 1, sistema%numeroTermalizacao
                call metropolis(T)
                !call wolff
                call superrelaxacao
            end do
           do i = 1, sistema%numeroPassosMC
                call metropolis(T)
                !call wolff
                call superrelaxacao
                call calcularSoma
            end do
            call cpu_time(tempoFinal)
            TempoCPU=tempoFinal-tempoInicial
            call calcularMedia
            t0=t0+dt
        end do
    end if

!fim do programa
CONTAINS
    !-----------------------------------------------------------------------------

    subroutine metropolis(T_)
        !Variaveis mudas 
        real, intent(in) :: T_
                
        !Variaveis Locais
        type(Spin) :: spinNovo
        type(Spin) ::  CampoEfetivo, aux_, aux2_
        integer :: i_, j_
        double precision :: deltaE 
        double precision :: probabilidade
        
        do i_=1, sistema%numeroSitios
            SpinNovo = direcao()
            deltaE = 0
            aux_ = campo(j_)
            CampoEfetivo%S(1) =aux_%S(1) +sistema%alfa*hr(j_)%h(1)
            CampoEfetivo%S(2) =aux_%S(2) + sistema%alfa*hr(j_)%h(2)
            do j_ = 1, 2
                deltaE = deltaE - (spinNovo%S(j_) - rede(i_)%S(j_))* CampoEfetivo%S(j_)
            end do
            probabilidade = 1.0d0/(1.0d0 + exp(deltaE/T_))
            call random_number(rando)
            if (rando <= probabilidade ) then
                !aceita a nova configuração
                forall  (j_ = 1:3)
                    quantidadesTermodinamicas%Magnetizacao(j_) = &
                        &quantidadesTermodinamicas%Magnetizacao(j_) + (spinNovo%S(j_) - rede(i_)%S(j_))
                end forall
                rede(i_) = SpinNovo
                quantidadesTermodinamicas%energia = quantidadesTermodinamicas%energia+deltaE
            end if

        end do 
    end Subroutine metropolis 
    !-----------------------------------------------------------------------------

    !    subroutine wolff ! errado
    !
    !        type :: Prato  !usado para construir a pilha
    !            double precision::prod
    !            integer::site
    !        end type Prato
    !
    !        !Variaveis Locais
    !        type(Prato), dimension(1 : sistema%numeroSitios) :: pilha
    !        type(Spin) :: direcaoSemente
    !        type(Spin) :: projecao
    !        type(Spin) :: projecaoVizinho
    !        type(Spin) :: campoEfetivo
    !        double precision :: produto
    !        double precision :: produtoSemente
    !        double precision :: produtoVizinho
    !        double precision :: produtoSitio
    !        double precision :: probabilidade
    !        double precision :: deltaE
    !        integer ::viz
    !        integer :: tamanhoCluster, semente, ponteiro, sitio, sitioVizinho
    !
    !        !inicia ponteiros do cluster
    !        ponteiro = 0 !diferente do programa velho porque eu incremento lá na frente
    !        tamanhoCluster = 0
    !
    !        !sorteia sitio semente
    !        do
    !            call random_number(Rando)
    !            semente = int(sistema%numeroSitios*Rando) + 1 !baseado no programa velho são sei se o gerador diferente pode prejudicar o resultado
    !            if (Rede(semente)%S2==1) exit
    !        end do
    !
    !        direcaoSemente = direcaoXY()
    !        produto = direcaoSemente%S(1)*rede(semente)%S(1) + direcaoSemente%S(2)*rede(semente)%S(2)
    !        produtoSemente = produto
    !
    !        !calcula projeção da semente na direção semente
    !        projecao%S(1) = direcaoSemente%S(1)*produto
    !        projecao%S(2) = direcaoSemente%S(2)*produto
    !
    !        !flipa semente
    !        rede(semente)%S(1) = rede(semente)%S(1) - 2.0d0*projecao%S(1)
    !        rede(semente)%S(2) = rede(semente)%S(2) - 2.0d0*projecao%S(2)
    !
    !        !Atualiza Energia e magnetização do sistema
    !        campoefetivo = campo(semente)
    !        deltaE = 2.0d0*(projecao%S(1)*campoefetivo%S(1) + projecao%S(2)*campoefetivo%S(2))
    !        quantidadesTermodinamicas%Energia = quantidadesTermodinamicas%Energia  + deltaE
    !        quantidadesTermodinamicas%Magnetizacao(1) = quantidadesTermodinamicas%Magnetizacao(1) -2.0d0*Projecao%S(1)
    !        quantidadesTermodinamicas%Magnetizacao(2) = quantidadesTermodinamicas%Magnetizacao(2) -2.0d0*Projecao%S(2)
    !
    !        !incrementa TamanhoCluster e cresce a pilha
    !        tamanhoCluster = tamanhoCluster + 1
    !        ponteiro = ponteiro + 1
    !        pilha(ponteiro)%site = semente
    !        pilha(ponteiro)%prod = produtoSemente
    !
    !        !repetição na pilha
    !        do while(ponteiro>0)
    !            !decrementa a pilha
    !            sitio = pilha(ponteiro)%site
    !            produtoSitio = pilha(ponteiro)%prod
    !            ponteiro = ponteiro - 1
    !            !repetição nos vizinhos
    !            do viz = 1, sistema%numeroCoordenacao
    !                sitioVizinho = Vizinhos(viz,sitio) !o primeiro marca o vizinho e segundo marca o sítio
    !                ! se o vizinho do sitio form magnetico continua
    !                if (rede(sitioVizinho)%S2==1) then
    !                    ProdutoVizinho = rede(sitioVizinho)%S(1)*DirecaoSemente%S(1) +rede(sitioVizinho)%S(2)*DirecaoSemente%S(2)
    !                    !se o vizinho está na mesma direção da semente
    !                    if (ProdutoVizinho*ProdutoSemente>0) then
    !                        probabilidade=1-exp(-2.0d0*ProdutoSitio*ProdutoVizinho/T) !!!ver necessidade do valor J
    !                        call random_number(Rando)
    !                        if (Rando<probabilidade) then
    !                            !aqui pode ser otimizado com forall
    !
    !                            !calcula projeção
    !                            ProjecaoVizinho%S(1)=produtoVizinho*DirecaoSemente%S(1)
    !                            ProjecaoVizinho%S(2)=produtoVizinho*DirecaoSemente%S(2)
    !
    !                            !flipa o spin
    !                            rede(sitioVizinho)%S(1) = rede(sitioVizinho)%S(1) - 2.0d0*projecaoVizinho%S(1)
    !                            rede(sitioVizinho)%S(2) = rede(sitioVizinho)%S(2) - 2.0d0*projecaoVizinho%S(2)
    !
    !                            !Atualiza Energia e magnetização do sistema
    !                            campoefetivo = campo(sitioVizinho)
    !                            deltaE = 2.0d0*(ProjecaoVizinho%S(1)*Campoefetivo%S(1) + ProjecaoVizinho%S(2)*campoefetivo%S(2))
    !                            quantidadesTermodinamicas%energia = quantidadesTermodinamicas%Energia  + deltaE
    !                            quantidadesTermodinamicas%magnetizacao(1) = &
    !                            &quantidadesTermodinamicas%magnetizacao(1) - 2.0d0*projecaoVizinho%S(1)
    !                            quantidadesTermodinamicas%magnetizacao(2) = &
    !                            &quantidadesTermodinamicas%magnetizacao(2) - 2.0d0*projecaoVizinho%S(2)
    !
    !                            !incrementa TamanhoCluster e cresce a pilha
    !                            tamanhoCluster = tamanhoCluster + 1
    !                            ponteiro = ponteiro + 1
    !                            pilha(ponteiro)%site = sitioVizinho
    !                            pilha(ponteiro)%prod = produtoVizinho
    !                       end if
    !                    end if
    !                end if
    !            end do  !repetição nos vizinhos
    !        end do !repetição na pilha
    !    end subroutine wolff
    !-----------------------------------------------------------------------------

    Subroutine superrelaxacao !ok
        !variaveis Locais
        integer :: i_, j_, k_
        type(spin) :: campoEfetivo, versorCampoEfetivo, novoSpin, projecao,aux_
        double precision :: moduloProjecao, moduloCampoEfetivo
            
   
        do i_ = 1 , sistema%numeroSitios
            call random_number(Rando)  
            j_ = int(sistema%numeroSitios*Rando )+ 1 !eu achava que deveria ser assim (NumeroSitios-1)*Rando+1  mas no programa aintigo é assim com está  não sei se há alguma diferença no numero aleatroio
            aux_ = campo(j_)
            CampoEfetivo%S(1) =aux_%S(1) +sistema%alfa*hr(j_)%h(1)
            CampoEfetivo%S(2) =aux_%S(2) + sistema%alfa*hr(j_)%h(2)

            moduloCampoEfetivo = sqrt(campoEfetivo%S(1)*campoEfetivo%S(1) + campoEfetivo%S(2)*campoEfetivo%S(2))
            versorCampoEfetivo%S(1) = campoEfetivo%S(1)/moduloCampoEfetivo
            versorCampoEfetivo%S(2) = campoEfetivo%S(2)/moduloCampoEfetivo
            moduloProjecao = versorCampoEfetivo%S(1)*Rede(J_)%S(1) + versorCampoEfetivo%S(2)*rede(J_)%S(2)
            forall (K_=1:2)
                projecao%S(k_) = moduloProjecao*versorCampoEfetivo%S(k_)
                novoSpin%S(k_) = - Rede(j_)%S(k_) + 2*projecao%S(k_)
                quantidadesTermodinamicas%magnetizacao(K_) = quantidadesTermodinamicas%magnetizacao(K_)&
                    & + novoSpin%S(K_) - Rede(J_)%S(K_)
                rede(j_)%S(k_) = novoSpin%S(K_)
            end forall

        end do
    End Subroutine Superrelaxacao
 
    !-----------------------------------------------------------------------------
 
    Subroutine calculaHelcidade
                     
        !Variaveis Locais
        integer :: i_ ,j_
        double precision :: helic1, helic2
         
   
        helic1 = 0.0d0
        helic2 = 0.0d0
        do i_ = 1, sistema%numeroSitios
            j_ = vizinhos(1,i_)
            helic1 = helic1 + rede(i_)%S(1)*rede(j_)%S(1) + Rede(i_)%S(2)*rede(j_)%S(2)
            helic2 = helic2 + rede(i_)%S(1)*rede(j_)%S(2) - Rede(i_)%S(2)*rede(j_)%S(1)
        end do
        quantidadesTermodinamicas%helicidade1 = helic1
        quantidadesTermodinamicas%helicidade2 = helic2
        quantidadesTermodinamicas%helicidade3 = helic2*helic2
    End Subroutine calculaHelcidade      
    !-----------------------------------------------------------------------------

    subroutine Inicializar !ok
        !Variaveis locais
        integer :: i_, j_
        type(Spin) :: campoefetivo,aux_,aux2_
            
        quantidadesTermodinamicas%magnetizacao=0 !magnetização instantanea
        quantidadesTermodinamicas%energia=0 !Energia instantanea
        do i_ = 1, sistema%numeroSitios
            aux_ = campo(i_)
            Campoefetivo%S(1) =aux_%S(1) + sistema%alfa*hr(i_)%h(1)
            Campoefetivo%S(2) =aux_%S(2) + sistema%alfa*hr(i_)%h(2)
            quantidadesTermodinamicas%energia = quantidadesTermodinamicas%energia - rede(i_)%S(1)*campoefetivo%S(1)
            quantidadesTermodinamicas%energia = quantidadesTermodinamicas%energia - rede(i_)%S(2)*campoefetivo%S(2)
        end do 
        forall  (j_=1:3)
            quantidadesTermodinamicas%magnetizacao(j_) = sum(rede(:)%S(j_))
        end forall
    end subroutine inicializar

    !-----------------------------------------------------------------------------
    subroutine marcarVizinhos !rede cubica !ok
         
        !Variaveis Locais
        integer,dimension(1:sistema%L) :: sucessor, antecessor

        integer :: i_, j_ , K_, L_
        
        L_=sistema%L    
     
        forall (i_ = 1:L_)
            sucessor(i_) = i_ + 1
            antecessor(i_) = i_ - 1
        end forall
        sucessor(L_) = 1
        antecessor(1) = L_
        forall (i_=1:L_, j_=1:L_, K_=1:L_)
            
            vizinhos(1,i_+(j_-1)*L_ + (K_-1)*L_*L_)=antecessor(i_) + (j_-1)*L_ + (K_-1)*L_*L_   
            vizinhos(2,i_+(j_-1)*L_ + (K_-1)*L_*L_)=sucessor(i_) + (j_-1)*L_ + (K_-1)*L_*L_
            
            vizinhos(3,i_+(j_-1)*L_ + (K_-1)*L_*L_)=i_ + (antecessor(j_)-1)*L_ + (K_-1)*L_*L_
            vizinhos(4,i_+(j_-1)*L_ + (K_-1)*L_*L_)=i_ + (sucessor(j_)-1)*L_ + (K_-1)*L_*L_
            
            vizinhos(5,i_+(j_-1)*L_ + (K_-1)*L_*L_)=i_ + (j_-1)*L_ + (antecessor(K_)-1)*L_*L_
            vizinhos(6,i_+(j_-1)*L_ + (K_-1)*L_*L_)=i_ + (j_-1)*L_ + (sucessor(K_)-1)*L_*L_
        end forall    

        
    end subroutine marcarVizinhos


    !   subroutine marcarsegVizinhos !rede cubica
    !        integer :: i_, K_
    !
    !        SegVizinhos=0
    !        do i_ = 1, sistema%numeroSitios
    !           k_ = vizinhos(1,i_)
    !           if (isolado(K_)==1) then
    !              segvizinhos(1,i_) = vizinhos(3,k_)
    !              segvizinhos(2,i_) = vizinhos(4,k_)
    !              segvizinhos(3,i_) = vizinhos(5,k_)
    !              segvizinhos(4,i_) = vizinhos(6,k_)
    !           end if
    !           k_ = vizinhos(2,i_)
    !           if (isolado(K_)==1) then
    !              segvizinhos(5,i_) = vizinhos(3,k_)
    !              segvizinhos(6,i_) = vizinhos(4,k_)
    !              segvizinhos(7,i_) = vizinhos(5,k_)
    !              segvizinhos(8,i_) = vizinhos(6,k_)
    !            end if
    !           k_ = vizinhos(3,i_)
    !           if (isolado(K_)==1) then
    !              segvizinhos(1,i_) = vizinhos(1,k_)
    !              segvizinhos(5,i_) = vizinhos(2,k_)
    !              segvizinhos(9,i_) = vizinhos(5,k_)
    !              segvizinhos(10,i_) = vizinhos(6,k_)
    !            end if
    !
    !            k_ = vizinhos(4,i_)
    !            if (isolado(K_)==1) then
    !              segvizinhos(2,i_) = vizinhos(1,k_)
    !              segvizinhos(6,i_) = vizinhos(2,k_)
    !              segvizinhos(11,i_) = vizinhos(5,k_)
    !              segvizinhos(12,i_) = vizinhos(6,k_)
    !            end if
    !
    !            k_ = vizinhos(5,i_)
    !            if (isolado(K_)==1) then
    !              segvizinhos(3,i_) = vizinhos(1,k_)
    !              segvizinhos(7,i_) = vizinhos(2,k_)
    !              segvizinhos(9,i_) = vizinhos(3,k_)
    !              segvizinhos(11,i_) = vizinhos(4,k_)
    !            end if
    !
    !            k_ = vizinhos(6,i_)
    !            if (isolado(K_)==1)then
    !              segvizinhos(4,i_) = vizinhos(1,k_)
    !              segvizinhos(8,i_) = vizinhos(2,k_)
    !              segvizinhos(10,i_) = vizinhos(3,k_)
    !              segvizinhos(12,i_) = vizinhos(4,k_)
    !            end if
    !
    !        end do
    !
    !
    !    end subroutine marcarsegVizinhos
    !-----------------------------------------------------------------------------
    subroutine iniciaVariaveis

        somaMag=0
        somaMag2=0
        somaMag4=0
        somaEne=0
        somaEne2=0
        somaEne4=0

    end subroutine iniciaVariaveis
    !-----------------------------------------------------------------------------
    Subroutine IniciaCampoAleatorio !ok
        !Variaveis Locais
        integer :: i_
        type(Spin) :: aux_
        do i_=1, sistema%numeroSitios
            aux_=direcaoXY
            hr(i_)%h(1)=aux_%S(1)
            hr(i_)%h(2)=aux_%S(2)
        end do

    end subroutine IniciaCampoAleatorio
    !-----------------------------------------------------------------------------

    subroutine lerDados(sistema_) !ok
    
        !Variaveis mudas 
        type(configuracao), intent(out) :: sistema_
                
        !Variaveis Locais
        character*20 :: lendo 
  
        read(1,*)  lendo 
        read(1,*)  sistema_%l
        write(*,*) lendo, sistema_%l
       
        read(1,*)  lendo 
        read(1,*)  sistema_%numeroCoordenacao
        write(*,*) lendo, sistema_%numeroCoordenacao
       
        read(1,*)  lendo 
        read(1,*)  sistema_%numeroTermalizacao
        write(*,*) lendo, sistema_%numeroTermalizacao
       
        read(1,*)  lendo 
        read(1,*)  sistema_%numeroPassosMC
        write(*,*) lendo, sistema_%numeroPassosMC
       
        read(1,*)  lendo 
        read(1,*)  sistema_%concentracao
        write(*,*) lendo, sistema_%concentracao
        
        read(1,*)  lendo 
        read(1,*)  sistema_%tInicial
        write(*,*) lendo, sistema_%tInicial
        
        read(1,*)  lendo 
        read(1,*)  sistema_%tFinal
        write(*,*) lendo, sistema_%tFinal
        
        read(1,*)  lendo 
        read(1,*)  sistema_%deltaT
        write(*,*) lendo, sistema_%deltaT
        sistema_%numeroSitios = sistema_%L*sistema_%L*sistema_%L
         
        read(1,*)  lendo 
        read(1,*)  sistema_%numeroAmostras
        write(*,*) lendo, sistema_%numeroAmostras
        
        read(1,*)  lendo 
        read(1,*)  sistema_%tsimulacao
        write(*,*) lendo, sistema_%tsimulacao
        
        read(1,*)  lendo 
        read(1,*)  sistema_%alfa !intensidade do campo aleatório
        write(*,*) lendo, sistema_%alfa
        write(*,*) "_____________________________"
        write(*,*) "_____________________________"
        sistema_%numeroSitios = sistema_%L*sistema_%L*sistema_%L
        sistema_%numeroIntervalosTeperatura = int((sistema_%tFinal -sistema_%tInicial)/sistema_%deltaT) +1
     
    
    end subroutine LerDados
    !!-----------------------------------------------------------------------------
    subroutine leargumentos(sistema_) !ok
          !Variaveis mudas
        type(configuracao), intent(out) :: sistema_




        character(len=32) :: arg
        Integer, Parameter :: wp = Selected_real_kind( 12, 70 )
        Real( wp ) :: a

        call get_command_argument(1, arg)
        read (arg,*) sistema_%l
        call get_command_argument(2, arg)
        read(arg,*) sistema_%alfa  !Intencidade de campo aleatório
        call get_command_argument(3, arg)
        read(arg,*) sistema_%tInicial
        call get_command_argument(4, arg)
        read(arg,*) sistema_%tFinal
        call get_command_argument(5, arg)
        read(arg,*) sistema_%deltaT
        call get_command_argument(6, arg)
        read(arg,*)sistema_%hinicial
        call get_command_argument(7, arg)
        read(arg,*) sistema_%hfinal
        call get_command_argument(8, arg)
        read(arg,*) sistema_%numeroTermalizacao
        call get_command_argument(9, arg)
        read(arg,*) sistema_%numeroPassosMC
        call get_command_argument(10, arg)
        read(arg,*) sistema_%histograma
        sistema_%numeroCoordenacao =6
    end subroutine leargumentos
    !-----------------------------------------------------------------------------
 
    subroutine gravarHistograma !ok
               
        !variaveis locais
        double precision :: energiaSitio, mx, my, mz
        integer :: numeroSitios_
        numeroSitios_ = sistema%numeroSitios
        energiaSitio = quantidadesTermodinamicas%Energia/NumeroSitios_
        mx = quantidadesTermodinamicas%Magnetizacao(1)/NumeroSitios_
        my = quantidadesTermodinamicas%Magnetizacao(2)/NumeroSitios_
        mz = quantidadesTermodinamicas%Magnetizacao(3)/NumeroSitios_
15      format (F18.7, 1x, F18.7, 1x, F18.7, 1x, F18.7)
        write(2,15) EnergiaSitio,Mx,My,Mz  
            
    end subroutine GravarHistograma
    !-----------------------------------------------------------------------------

    function direcaoXY() !ok
        !variaveis locais
        type(spin) :: direcaoXY
        double precision :: theta
        double precision :: PI=3.1415926535897932384626433832795028841971693993751 
        
        call random_number(rando)
        theta = PI*rando
        direcaoXY%S(1) = cos(Theta)
        direcaoXY%S(2) = sin(Theta)
        direcaoXY%S(3) = 0

    end function direcaoXY 
    !-----------------------------------------------------------------------------
   
    function campo(i_)
        
        !Variaveis mudas 
        type(Spin) :: campo 
        integer, intent(in) :: i_
        
        !variaveis locais
        integer :: j_,K_
        
        campo%S=(/0,0,0/)
        do j_ = 1, sistema%numeroCoordenacao
            k_ = vizinhos(j_,i_)
            campo%S(1) = campo%S(1) + rede(K_)%S(1) 
            campo%S(2) = campo%S(2) + rede(K_)%S(2)
        end do

    
    end function campo    
    !-----------------------------------------------------------------------------

    function direcao()!Marsaglia(rand) !ok
        
        !Variaveis mudas 
        type(spin) :: direcao
       
        !variaveis locais
        double precision :: auxiliar1, auxiliar2, auxiliar3, auxiliar4  
        
        do 
            call random_number(rando)
            auxiliar1 = 1- 2*rando
            call random_number(Rando)
            auxiliar2= 1- 2*Rando
            auxiliar3 = auxiliar1*auxiliar1 + auxiliar2*auxiliar2
            if (auxiliar3 <= 1) exit
        end do
        auxiliar4 = sqrt(1 - auxiliar3)
        direcao%S(1) = 2*auxiliar1*auxiliar4
        direcao%S(2) = 2*auxiliar2*auxiliar4
        direcao%S(3) = 1 - 2*auxiliar3

    end function direcao !Marsaglia
 !-----------------------------------------------------------------------------
        subroutine calcularMedia
        integer ::numeroSitios,L, MCc
        L=sistema%l
        MCc=sistema%numeroPassosMC
        numeroSitios = 2*L**3
        mediaMag = somaMag/MCc
        mediaMag2 = somaMag2/MCc
        mediaMag4 = somaMag4/MCc
        mediaEne=somaEne/MCc
        mediaEne2=somaEne2/MCc
        mediaEne4=somaEne4/MCc
        calor= (MediaEne2-mediaEne*mediaEne)/NumeroSitios/t0/t0
        susceptibilidade= (mediaMag2 - MediaMag*MediaMag)*NumeroSitios/t0
        cumuM=1-mediaMag4/(3*mediaMag2*mediaMag2)  ! rever essa equaÃ§Ã£o soma ou media
        write(*,*)t0, susceptibilidade, calor, cumuM, mediaMag, mediaEne, TempoCPU
    end subroutine calcularMedia
    !-----------------------------------------------------------------------------
      subroutine calcularSoma
        double precision ::mag, mag2, mag4, ene,ene2, ene4,L
        L=sistema%l
        mag=quantidadesTermodinamicas%Magnetizacao(1)/(2*L**3)
        mag2=mag*mag
        mag4=mag2*mag2
        ene =  quantidadesTermodinamicas%Energia
        ene2 = ene*ene
        ene4 = ene2*ene2
        somaMag=somaMag+mag
        somaMag2=somaMag2+mag2
        somaMag4=somaMag4+mag4
        somaEne=somaEne+ene
        somaEne2=somaEne2+ene2
        somaEne4=somaEne4+ene4

    end subroutine calcularSoma
    !--------------------

    !-----------------------------------------------------------------------------
  
    SUBROUTINE init_random_seed()
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
    END SUBROUTINE
end program 
