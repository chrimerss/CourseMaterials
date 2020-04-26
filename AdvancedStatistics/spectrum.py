class PowerSpectrum(object):
    '''
    Power spectrum analysis with three techniques and significance test
    
    Examples:
    --------------------------
    #with normal condition
    ps= PowerSpectrum(y,dt=1,detrend=True,window=False)
    x, ps_curve, rn_curve, sig_curve= ps.normal()
    
    #with running mean (smooth_window=5)
    ps= PowerSpectrum(y,dt=1,detrend=True,window=False)
    x, ps_curve, rn_curve, sig_curve= ps.running_mean(window=5)
    
    #with window (hann window) tapering
    ps= PowerSpectrum(y,dt=1,detrend=True,window=True, subset=5)
    x, ps_curve, rn_curve, sig_curve= ps.window(fw=1.2)
    '''
    
    def __init__(self, y, dt=1, detrend=True, window=False, **window_params):
        '''
        Args:
        ---------------------
        :y - np.ndarray; geophysical values (anomaly)
        :dt - int; time step;
        :detrend - bool; whether to detrend your data;
        :window - bool; whether to adopt taper and subdivision
        :window_params: dict {'subset': int,
                              ...
                              }
                              
        '''
        if detrend:
            y,_ ,_ = self._detrend(y)
        num= len(y)
        self.N = int(np.floor(num/2.)*2)
        self.y = y[:self.N]
        self.nyquist= int(self.N/2)
        self.T= dt * self.N
        i = np.arange(1,self.N+1)
        if window:
            self.subset= window_params.get('subset')
            segments= self.subset*2-1
            num = len(y)/self.subset #The length of the time series over which to conduct power spectrum analysis
            self.N = int(np.floor(num/2.)*2)
            self.nyquist= int(self.N/2)
            self.T= dt * self.N
            i = np.arange(1,self.N+1)
            A = np.ones((self.nyquist+1,segments))*np.nan
            B = np.ones((self.nyquist+1,segments))*np.nan
            self.rho= []
            w= signal.windows.hann(self.N)
            #--------compute each subset-----------#
            for n in range(self.subset):
                ysub = y[n*self.N:(n+1)*self.N]
                self.rho.append(stats.pearsonr(ysub[1:], ysub[:-1])[0])
                ysub*= w

                for k in range(1, self.nyquist): #resolve the lowest frequency
                    A[k,n] = 2./self.N*(ysub*np.cos(2*np.pi*k*i*dt/self.T)).sum()
                    B[k,n] = 2./self.N*(ysub*np.sin(2*np.pi*k*i*dt/self.T)).sum()
                k= self.nyquist
                A[k,n] = 1/self.N*(ysub*np.cos(np.pi*self.N*i*dt/self.T)).sum()
                B[k,n]= 0 
            #---------compute overlaping----------#
            start= self.N//2
            for n in range(self.subset, segments):
                ysub= y[start:start+self.N]
#                 print(start,start+self.N)
                self.rho.append(stats.pearsonr(ysub[1:], ysub[:-1])[0])
                ysub*= w
                for k in range(1, self.nyquist): #resolve the lowest frequency
                    A[k,n] = 2./self.N*(ysub*np.cos(2*np.pi*k*i*dt/self.T)).sum()
                    B[k,n] = 2./self.N*(ysub*np.sin(2*np.pi*k*i*dt/self.T)).sum()
                k= self.nyquist
                A[k,n] = 1/self.N*(ysub*np.cos(np.pi*self.N*i*dt/self.T)).sum()
                B[k,n]= 0 
                start+=self.N
            C2= (A**2+ B**2)[1:]

            self.C2 = np.nanmean(C2.copy(), axis=1)
            
        else:
            A = np.ones((self.nyquist+1,))*np.nan
            B = np.ones((self.nyquist+1,))*np.nan
            for k in range(1,self.nyquist):
                A[k] = 2./self.N*(self.y*np.cos(2*np.pi*k*i*dt/self.T)).sum()
                B[k] = 2./self.N*(self.y*np.sin(2*np.pi*k*i*dt/self.T)).sum()
            #-------------------------------------------------------------------
            # Calculate A_k and B_k at k = nyquist.
            #-------------------------------------------------------------------
            k = self.nyquist
            A[k] = 1/self.N*(self.y*np.cos(np.pi*self.N*i*dt/self.T)).sum()
            B[k] = 0

            #-------------------------------------------------------------------
            # Calculate C_k**2 - i.e., the total magnitude.
            #-------------------------------------------------------------------
            C2 = A**2 + B**2
            self.C2 = C2[1:].copy() # Keep only the k = 1, 2, 3,... values. Skip the k=0 wave.
            self.rho= stats.pearsonr(self.y[1:], self.y[:-1])[0]
    
    def normal(self, priori=False):
        '''
        Return the scaled results for angular frequency, power specturm density, red noise power specturm
        and significance curve for plot
        '''
        C2= self.C2
        dof=2
        rho= self.rho
        return self.sig_test(C2, dof, rho, priori)
    
    def running_mean(self,window, priori=False):
        '''
        Return the scaled results for angular frequency, power specturm density, red noise power specturm
        and significance curve for plot
        '''        
        dof=2
        RUNMEAN = window
        C2 = np.convolve(self.C2, np.ones((RUNMEAN,))/RUNMEAN)[(RUNMEAN-1):] #Fast way of doing a running mean.
        dof*=RUNMEAN # Applying a running mean increases the degrees of freedom for significance
                         # testing.
        rho= np.nanmean(self.rho)
        
        return self.sig_test(C2, dof, rho, priori)
    
    def window(self,fw=1.2, priori=False):
        '''
        Return the scaled results for angular frequency, power specturm density, red noise power specturm
        and significance curve for plot
        '''        
        C2= self.C2
        fw= fw
        dof= np.ceil(2*self.subset*fw)
        rho= np.nanmean(self.rho)
        
        return self.sig_test(C2, dof, rho, priori)
    
    
    def sig_test(self, C2, dof, rho, priori=False):
        #significance test
        k = np.arange(1,self.nyquist+1) # Wavenumbers
        omega = 2*np.pi*k/self.N # Angular frequency

        sigma_e2 = (1-rho**2)*np.var(self.y)

        rnPower = 4*sigma_e2**2/self.N/(1+rho**2-2*rho*np.cos(omega))
        delta = omega[1]-omega[0]
        scaledPower = C2/np.nansum(delta*C2)
        scaledRN = rnPower/np.nansum(delta*rnPower)
        
        alphaStar = 0.05
        if priori: alpha=alphaStar
        else: alpha = alphaStar/self.nyquist

        sigChi2 = stats.chi2.isf(alpha,dof) # Chi-squared value for the nominal alpha value for your 95% test.

        sigCurve = scaledRN / dof *sigChi2 
        
        return omega, scaledPower, scaledRN, sigCurve
    
    def _detrend(self,y):
        '''
        linear regression to detrend

        Args:
        -------------------
        :y - numpy.ndarray; target to detrend

        Returns:
        -------------------
        :detrend - numpy.ndarray; 
        '''

        t= np.arange(len(y))
        E= np.ones((t.size,2))*np.nan
        E[:,0] = t
        E[:,1] = 1
        xhat = np.linalg.inv(np.transpose(E).dot(E)).dot(E.T).dot(y)
        trend = E.dot(xhat)
        detrend= y- trend

        return detrend, xhat, trend
