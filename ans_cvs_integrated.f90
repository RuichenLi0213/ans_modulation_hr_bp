! ================================================================
! ans_cvs_integrated.f90
!
! Integrated reduced ANS-CVS model. Combines:
!   - Neuro_HR formula extracted from ans_modulation.f90
!     (which itself follows the structure of Neuro_HR in cvs.f90,
!      with the physiologically corrected sigmoid sign convention)
!   - Neuro_BP formula from bp_modulation_reflected.f90
!     (deviation form of the cvs.f90 Neuro_BP)
!   - Half-sine cardiac inflow + delayed reflected inflow with
!     mass-preservation rescaling, from bp_modulation_reflected.f90
!   - A simplified Disturbance() subroutine that acts only on the
!     vascular-bed resistance, mimicking the role of Disturbance in
!     cvs.f90 (which set SVB_R_new = 2*R0 during a window).
!
! What is preserved verbatim (symbolic structure):
!   ns_hr = 1/(1+(hr/miu_hr)**niu_s_hr)
!   np_hr = 1/(1+(hr/miu_hr)**niu_p_hr)
!   EneuHR_hr = AneuHR_hr*ns_hr - BneuHR_hr*np_hr + CneuHR_hr
!   RK4 update of HR
!   ns_p = 1/(1+(P_afferent/miu_p)**niu_s_p)
!   np_p = 1/(1+(P_afferent/miu_p)**niu_p_p)
!   R_eq = R0_base + AneuR*(ns_p-0.5) - BneuR*(np_p-0.5)
!   First-order Euler relaxation of R_vb
!   Q_in = Q_forward + Q_reflected(t - tau_reflect)
!   Mass-preserving Q_peak = SV * pi / (2 * t_systole * (1 + k_reflect))
!
! Interface adaptations (not changes in meaning):
!   tneuhr_hr is documented as time constant in SECONDS, matching
!   the unit convention used for TneuR and tau_aff in the BP model.
!   The afferent BP signal is the low-pass-filtered P_art (Pfa-like).
!   Disturbance acts on a single R_vb instead of 11 vascular beds.
!
! Units:
!   pressure mmHg, time s, flow mL/s, volume mL,
!   compliance mL/mmHg, resistance mmHg*s/mL, HR beats/s.
! ================================================================
module ans_cvs_globals
  implicit none

  ! ---------------- Simulation control ----------------
  integer      :: nstep, nlast, nprint
  integer      :: n_onset           ! step at which both ANS modules engage
  integer      :: n_dist_onset      ! disturbance start step
  integer      :: n_dist_end        ! disturbance end step
  real(kind=8) :: dt
  real(kind=8) :: time

  ! ---------------- HR state (from ans_modulation.f90) ----------------
  real(kind=8) :: tduration         ! current cardiac cycle (s)
  real(kind=8) :: tduration_new
  real(kind=8) :: hr                ! current HR (beats/s)
  real(kind=8) :: hr_new1, hr_new2
  real(kind=8) :: hr0               ! initial HR

  ! ---------------- HR neuromodulation parameters ----------------
  real(kind=8) :: miu_hr            ! target HR (beats/s)
  real(kind=8) :: niu_s_hr, niu_p_hr
  real(kind=8) :: ns_hr, np_hr
  real(kind=8) :: AneuHR_hr, BneuHR_hr, CneuHR_hr
  real(kind=8) :: TneuHR_hr         ! time constant (s)
  real(kind=8) :: EneuHR_hr         ! equilibrium HR (after ans_gain blend)
  real(kind=8) :: EneuHR_raw        ! raw equilibrium HR from the verbatim formula

  ! ---------------- Cardiac inflow + reflected wave ----------------
  real(kind=8) :: t_systole         ! systolic ejection duration (s)
  real(kind=8) :: SV                ! stroke volume (mL, fixed)
  real(kind=8) :: tcr, tcr_prev     ! time within current cycle
  ! Fix: incremental cardiac-phase accumulator. The phase is advanced by
  ! dt each step and wrapped against the CURRENT tduration, instead of
  ! being rebuilt from absolute time as time - tduration*int(time/tduration).
  ! Rebuilding from absolute time makes tcr jump discontinuously whenever
  ! tduration changes (HR modulation), which relocates the ejection phase
  ! and injects a spurious inflow pulse -> the BP spike at ANS onset.
  logical      :: cycle_wrapped     ! .true. on the step a cycle restarts
  real(kind=8) :: Q_peak
  real(kind=8) :: Q_forward, Q_reflected, Q_in

  real(kind=8) :: k_reflect_base    ! reflection strength
  real(kind=8) :: k_reflect_eff
  real(kind=8) :: tau_reflect       ! reflection delay (s)
  integer      :: ndelay            ! delay in steps
  integer      :: nbuf
  real(kind=8), allocatable :: Q_fwd_buf(:)
  logical      :: mass_preserve     ! .true. keeps cycle inflow = SV
  logical      :: couple_kreflect_R ! .true. couples k_reflect to R_vb

  ! ---------------- Vascular state and parameters ----------------
  real(kind=8) :: P_art, P_ven
  real(kind=8) :: Q_out, Q_return
  real(kind=8) :: C_art, C_ven
  real(kind=8) :: R_vb              ! current vascular-bed resistance
  real(kind=8) :: R0_base           ! unperturbed baseline R
  real(kind=8) :: R_ven, P_RA

  ! ---------------- BP neuromodulation ----------------
  real(kind=8) :: P_afferent        ! low-pass-filtered P_art (Pfa-like)
  real(kind=8) :: tau_aff           ! afferent filter time const (s)
  real(kind=8) :: miu_p             ! target arterial pressure (mmHg)
  real(kind=8) :: niu_s_p, niu_p_p
  real(kind=8) :: ns_p, np_p
  real(kind=8) :: AneuR, BneuR, TneuR
  real(kind=8) :: R_eq              ! equilibrium R (after ans_gain blend)
  real(kind=8) :: R_eq_raw          ! raw equilibrium R from the verbatim formula

  ! ---------------- Smooth ANS onset (fix for BP spike) ----------------
  ! ans_gain ramps 0 -> 1 over T_ramp seconds after t_onset. It multiplies
  ! ONLY the ANS correction (deviation) terms, leaving the verbatim
  ! Neuro_HR / Neuro_BP formulas unchanged.
  real(kind=8) :: ans_gain
  real(kind=8) :: T_ramp            ! ANS activation ramp duration (s)
  real(kind=8) :: t_onset           ! onset time in seconds (= n_onset*dt)

  ! ---------------- Disturbance ----------------
  real(kind=8) :: R_dist_factor     ! multiplier applied to R0_base
                                    ! during the disturbance window
  logical      :: disturbance_active

  ! ---------------- Cycle summary ----------------
  integer      :: cycle_index
  real(kind=8) :: P_art_max_cycle, P_art_min_cycle
  real(kind=8) :: P_art_sum_cycle, P_art_count_cycle
  real(kind=8) :: SBP, DBP, MAP_cycle, pulse_pressure
end module ans_cvs_globals


! ================================================================
! Main driver
! ================================================================
program ans_cvs_integrated
  use ans_cvs_globals
  implicit none

  call parameterset()

  ! Build the delay buffer based on tau_reflect set in parameterset()
  ndelay = max(1, int(tau_reflect / dt))
  nbuf   = ndelay + 10
  allocate(Q_fwd_buf(nbuf))
  Q_fwd_buf = 0.0d0

  do nstep = 1, nlast
    time = dt * dble(nstep)
    call neuro_hr()
    call cardiac_inflow()
    call neuro_bp()
    call disturbance()
    call update_cvs()
    call update_summary()
    call output_data()
  end do

  deallocate(Q_fwd_buf)
end program ans_cvs_integrated


! ================================================================
! HR neuromodulation. Symbolic form preserved exactly from
! ans_modulation.f90 (which is the simplified form of Neuro_HR
! from cvs.f90, with corrected sigmoid signs).
!
! Interface adaptation: TneuHR_hr is now a time constant in seconds.
! ================================================================
subroutine neuro_hr()
  use ans_cvs_globals
  implicit none
  real(kind=8) :: k1, k2, k3, k4
  real(kind=8), parameter :: pi = 3.14159265358979d0

  ! Setpoint, sigmoid slopes (verbatim from ans_modulation.f90)
  miu_hr   =  1.25d0    ! 75 bpm
  niu_s_hr =  7.0d0
  niu_p_hr = -7.0d0

  ! Gains and baseline offset (verbatim from ans_modulation.f90)
  AneuHR_hr = 1.0d0
  BneuHR_hr = 1.0d0
  CneuHR_hr = 1.333d0
  TneuHR_hr = 5.0d0     ! seconds (cvs.f90 used 100000 steps = 5 s)

  ! Fix: smooth ANS onset to avoid non-physiological BP spike.
  ! ans_gain ramps 0 -> 1 with a raised-cosine over T_ramp seconds.
  ! It is applied ONLY to the deviation of EneuHR_raw from its own
  ! baseline CneuHR_hr; the verbatim sigmoid formula is unchanged.
  if (nstep < n_onset) then
    ans_gain = 0.0d0
  else
    ans_gain = 0.5d0 * (1.0d0 - cos(pi * min(1.0d0, &
               (time - t_onset) / T_ramp)))
  end if

  ! Sigmoid activities computed EVERY step from the current HR, so they
  ! are continuous through onset (no hard switch from a fixed 0.5).
  ns_hr = 1.0d0 / (1.0d0 + (hr_new1 / miu_hr) ** niu_s_hr)
  np_hr = 1.0d0 / (1.0d0 + (hr_new1 / miu_hr) ** niu_p_hr)

  ! Verbatim equilibrium-HR formula from ans_modulation.f90
  EneuHR_raw = AneuHR_hr * ns_hr - BneuHR_hr * np_hr + CneuHR_hr

  ! Fix: blend toward the raw equilibrium through ans_gain so that at
  ! onset (ans_gain = 0) the target equals the formula baseline CneuHR_hr,
  ! giving zero initial driving deviation and hence no HR jump.
  EneuHR_hr = CneuHR_hr + ans_gain * (EneuHR_raw - CneuHR_hr)

  ! Continuous RK4 update of HR toward EneuHR_hr every step (no branch).
  ! Before onset ans_gain = 0 so EneuHR_hr = CneuHR_hr; with hr starting
  ! at 100 bpm this would slowly pull HR toward the baseline, so to keep
  ! HR exactly fixed during warm-up we hold it until onset.
  if (nstep < n_onset) then
    hr_new2 = hr_new1
  else
    k1 = dt * (-hr_new1                + EneuHR_hr) / TneuHR_hr
    k2 = dt * (-(hr_new1 + 0.5d0 * k1) + EneuHR_hr) / TneuHR_hr
    k3 = dt * (-(hr_new1 + 0.5d0 * k2) + EneuHR_hr) / TneuHR_hr
    k4 = dt * (-(hr_new1 +         k3) + EneuHR_hr) / TneuHR_hr
    hr_new2 = hr_new1 + (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) / 6.0d0
  end if
  tduration_new = 1.0d0 / hr_new2

  ! Commit the new HR for downstream use this same step
  tduration = tduration_new
  hr_new1   = 1.0d0 / tduration
  hr        = hr_new1
end subroutine neuro_hr


! ================================================================
! Cardiac inflow with phenomenological wave reflection.
! Half-sine forward ejection of stroke volume SV during systole,
! plus a delayed, scaled copy from a circular buffer.
! Uses the CURRENT dynamic HR set by neuro_hr().
! ================================================================
subroutine cardiac_inflow()
  use ans_cvs_globals
  implicit none
  real(kind=8), parameter :: pi = 3.14159265358979d0
  real(kind=8) :: k_scale
  integer :: i_write, i_read

  ! Fix: advance the cardiac phase incrementally and wrap against the
  ! CURRENT tduration. This keeps tcr (and therefore Q_forward) continuous
  ! when HR changes, eliminating the phase jump that caused the BP spike.
  tcr = tcr + dt
  if (tcr >= tduration) then
    tcr = tcr - tduration       ! restart cycle, carrying the remainder
    cycle_wrapped = .true.
  else
    cycle_wrapped = .false.
  end if

  ! Optional ANS->reflection coupling
  if (couple_kreflect_R) then
    k_reflect_eff = k_reflect_base * (R_vb / R0_base)
  else
    k_reflect_eff = k_reflect_base
  end if

  ! Mass-preservation rescaling
  if (mass_preserve) then
    k_scale = 1.0d0 + k_reflect_eff
  else
    k_scale = 1.0d0
  end if
  Q_peak = SV * pi / (2.0d0 * t_systole * k_scale)

  ! Forward ejection
  if (tcr < t_systole) then
    Q_forward = Q_peak * sin(pi * tcr / t_systole)
  else
    Q_forward = 0.0d0
  end if

  ! Write into the delay buffer
  i_write = mod(nstep - 1, nbuf) + 1
  Q_fwd_buf(i_write) = Q_forward

  ! Read the delayed value
  if (nstep > ndelay) then
    i_read = mod(nstep - 1 - ndelay, nbuf) + 1
    Q_reflected = k_reflect_eff * Q_fwd_buf(i_read)
  else
    Q_reflected = 0.0d0
  end if

  Q_in = Q_forward + Q_reflected
end subroutine cardiac_inflow


! ================================================================
! BP neuromodulation. Symbolic form preserved from
! bp_modulation_reflected.f90 (deviation form of cvs.f90 Neuro_BP).
! ================================================================
subroutine neuro_bp()
  use ans_cvs_globals
  implicit none
  real(kind=8) :: dP_aff, dR

  ! Low-pass filter P_art to obtain afferent (MAP-like) signal
  dP_aff     = (P_art - P_afferent) / tau_aff
  P_afferent = P_afferent + dt * dP_aff

  ! Setpoint and sigmoid slopes (verbatim from bp_modulation_reflected.f90)
  miu_p   = 95.0d0
  niu_s_p =  7.0d0
  niu_p_p = -7.0d0

  ! Gains and time constant (verbatim from bp_modulation_reflected.f90)
  AneuR = 0.6d0
  BneuR = 0.0d0
  TneuR = 10.0d0    ! seconds

  ! Fix: compute the sigmoid activities EVERY step (continuous through
  ! onset) rather than hard-switching from a fixed 0.5. The verbatim
  ! deviation formula is preserved; ans_gain blends the raw equilibrium
  ! toward R0_base so that at onset (ans_gain = 0) R_eq = R0_base exactly,
  ! giving zero driving deviation and no R_vb / BP jump.
  ns_p = 1.0d0 / (1.0d0 + (P_afferent / miu_p) ** niu_s_p)
  np_p = 1.0d0 / (1.0d0 + (P_afferent / miu_p) ** niu_p_p)

  ! Verbatim equilibrium-R formula from bp_modulation_reflected.f90
  R_eq_raw = R0_base + AneuR * (ns_p - 0.5d0) - BneuR * (np_p - 0.5d0)

  ! Smooth-onset blend (ans_gain already set in neuro_hr this step)
  R_eq = R0_base + ans_gain * (R_eq_raw - R0_base)

  ! First-order continuous relaxation of R_vb toward R_eq every step.
  ! During warm-up ans_gain = 0 so R_eq = R0_base; we still hold R_vb so
  ! that disturbance() retains full control of the warm-up resistance.
  if (nstep < n_onset) then
    ! R_vb left to disturbance()/initial value; no ANS relaxation yet.
    continue
  else
    dR   = (R_eq - R_vb) / TneuR
    R_vb = R_vb + dt * dR
  end if
end subroutine neuro_bp


! ================================================================
! Disturbance. Simplified version of cvs.f90 Subroutine Disturbance.
! Acts ONLY on the vascular-bed resistance: during the window
!   n_dist_onset <= nstep <= n_dist_end
! the resistance is forced to R_dist_factor * R0_base, mimicking the
! original "SVB_R_new(i) = 2*R0(i)" perturbation. Outside the window,
! resistance is left to the current state (either R0_base during the
! pre-ANS phase, or whatever Neuro_BP has set after onset).
! Disturbance does NOT modify HR, SV, compliance, or any other state.
! ================================================================
subroutine disturbance()
  use ans_cvs_globals
  implicit none

  if (nstep >= n_dist_onset .and. nstep <= n_dist_end) then
    R_vb = R_dist_factor * R0_base
    disturbance_active = .true.
  else
    disturbance_active = .false.
  end if
end subroutine disturbance


! ================================================================
! Closed-loop CVS update by explicit Euler.
!   Q_out    = (P_art - P_ven) / R_vb
!   Q_return = (P_ven - P_RA)  / R_ven       (clamped >= 0)
!   dP_art/dt = (Q_in - Q_out)    / C_art
!   dP_ven/dt = (Q_out - Q_return) / C_ven
! ================================================================
subroutine update_cvs()
  use ans_cvs_globals
  implicit none
  real(kind=8) :: dPa, dPv

  Q_out    = (P_art - P_ven) / R_vb
  Q_return = (P_ven - P_RA)  / R_ven
  if (Q_return < 0.0d0) Q_return = 0.0d0

  dPa = (Q_in  - Q_out)    / C_art
  dPv = (Q_out - Q_return) / C_ven
  P_art = P_art + dt * dPa
  P_ven = P_ven + dt * dPv
end subroutine update_cvs


! ================================================================
! Per-cycle SBP / DBP / MAP / pulse-pressure tracking. A cycle ends
! when tcr wraps back near zero (tcr < tcr_prev).
! ================================================================
subroutine update_summary()
  use ans_cvs_globals
  implicit none

  if (P_art > P_art_max_cycle) P_art_max_cycle = P_art
  if (P_art < P_art_min_cycle) P_art_min_cycle = P_art
  P_art_sum_cycle   = P_art_sum_cycle   + P_art
  P_art_count_cycle = P_art_count_cycle + 1.0d0

  ! Detect a genuine cycle boundary using the incremental phase wrap flag
  ! set in cardiac_inflow(). This is robust to HR changes (no false
  ! triggers from absolute-time phase reconstruction).
  if (cycle_wrapped) then
    SBP = P_art_max_cycle
    DBP = P_art_min_cycle
    if (P_art_count_cycle > 0.0d0) then
      MAP_cycle = P_art_sum_cycle / P_art_count_cycle
    else
      MAP_cycle = P_art
    end if
    pulse_pressure = SBP - DBP
    cycle_index = cycle_index + 1

    write(13, '(i8, 7f16.6)') cycle_index, time, hr, SBP, DBP, &
                              MAP_cycle, pulse_pressure, R_vb

    P_art_max_cycle   = P_art
    P_art_min_cycle   = P_art
    P_art_sum_cycle   = 0.0d0
    P_art_count_cycle = 0.0d0
  end if
  tcr_prev = tcr
end subroutine update_summary


! ================================================================
! Output files (headers written at nstep == 1):
!   ANS_CVS_waveform.dat
!     time, HR, P_art, P_ven, Q_forward, Q_reflected,
!     Q_in_total, Q_out, Q_return, R_vb
!   ANS_activity.dat
!     time, HR_afferent, BP_afferent, ns_hr, np_hr,
!     ns_bp, np_bp, HR_eq, R_eq, HR, R_vb
!   ANS_CVS_summary.dat   (written from update_summary())
!     cycle, time, HR, SBP, DBP, MAP, pulse_pressure, R_vb
! ================================================================
subroutine output_data()
  use ans_cvs_globals
  implicit none

  if (nstep == 1) then
    open(11, file='ANS_CVS_waveform.dat', status='replace')
    open(12, file='ANS_activity.dat',     status='replace')
    open(13, file='ANS_CVS_summary.dat',  status='replace')

    write(11,'(A)') &
      '# time  HR  P_art  P_ven  Q_forward  Q_reflected  Q_in_total  Q_out  Q_return  R_vb'
    write(12,'(A)') &
      '# time  HR_afferent  BP_afferent  ns_hr  np_hr  ns_bp  np_bp  HR_eq  R_eq  HR  R_vb  ans_gain'
    write(13,'(A)') &
      '# cycle  time  HR  SBP  DBP  MAP  pulse_pressure  R_vb'
  end if

  if (mod(nstep, nprint) == 0 .or. nstep == 1) then
    write(11,'(10f16.6)') time, hr, P_art, P_ven, Q_forward, &
                          Q_reflected, Q_in, Q_out, Q_return, R_vb
    write(12,'(12f16.6)') time, hr, P_afferent, ns_hr, np_hr, &
                          ns_p, np_p, EneuHR_hr, R_eq, hr, R_vb, ans_gain
  end if
end subroutine output_data


! ================================================================
! Centralized parameter and initial-condition setup.
! ================================================================
subroutine parameterset()
  use ans_cvs_globals
  implicit none
  real(kind=8) :: dx, crno, c_wave

  ! ---------------- Simulation control ----------------
  nprint        = 100
  nlast         = 1800000     ! ~90 s total
  n_onset       =  400000     ! ANS modules engage at 20 s
  n_dist_onset  =  100000     ! disturbance window:  5 s
  n_dist_end    =  300000     !                     to 15 s (before ANS onset)

  ! Time step (same convention as ans_modulation.f90 / bp_modulation_reflected.f90)
  dx     = 0.001d0
  crno   = 0.5d0
  c_wave = 10.0d0
  dt     = dx * crno / c_wave        ! = 5.0e-5 s

  ! ---------------- HR initial state ----------------
  tduration     = 0.6d0              ! 100 bpm
  tduration_new = tduration
  hr            = 1.0d0 / tduration
  hr0           = hr
  hr_new1       = hr
  hr_new2       = hr

  ! ---------------- Cardiac inflow ----------------
  t_systole = 0.30d0                 ! systolic ejection (s)
  SV        = 70.0d0                 ! stroke volume (mL, FIXED)
  tcr       = 0.0d0
  tcr_prev  = 0.0d0
  cycle_wrapped = .false.

  ! ---------------- Reflected wave ----------------
  k_reflect_base    = 0.30d0
  tau_reflect       = 0.25d0
  mass_preserve     = .true.
  couple_kreflect_R = .false.

  ! ---------------- Vascular ----------------
  C_art  = 1.5d0
  C_ven  = 100.0d0
  R0_base = 1.0d0                    ! unperturbed baseline R
  R_vb   = R0_base                   ! initial R = baseline
  R_ven  = 0.04d0
  P_RA   = 2.0d0

  ! ---------------- Disturbance ----------------
  R_dist_factor = 1.5d0              ! mild disturbance; set to 2.0 to mimic
                                     ! cvs.f90 exactly, or 1.0 to disable
  disturbance_active = .false.

  ! ---------------- Initial pressures ----------------
  P_art      = 80.0d0
  P_ven      =  5.0d0
  P_afferent = P_art
  tau_aff    =  2.0d0

  ! ---------------- ANS default neutral values ----------------
  ns_hr     = 0.5d0;  np_hr = 0.5d0
  ns_p      = 0.5d0;  np_p  = 0.5d0
  EneuHR_hr = hr0
  EneuHR_raw = hr0
  R_eq      = R0_base
  R_eq_raw  = R0_base

  ! ---------------- Smooth ANS onset (fix for BP spike) ----------------
  ! ans_gain ramps 0 -> 1 over T_ramp seconds after t_onset. Applied only
  ! to the ANS deviation terms; the raw Neuro_HR / Neuro_BP formulas are
  ! unchanged. T_ramp ~ 3-5 s gives a smooth, physiological engagement.
  t_onset  = dt * dble(n_onset)      ! = 20 s with the defaults
  T_ramp   = 4.0d0                   ! ANS activation ramp duration (s)
  ans_gain = 0.0d0

  ! ---------------- Cycle tracking ----------------
  cycle_index       = 0
  P_art_max_cycle   = P_art
  P_art_min_cycle   = P_art
  P_art_sum_cycle   = 0.0d0
  P_art_count_cycle = 0.0d0
  SBP = P_art;  DBP = P_art;  MAP_cycle = P_art
  pulse_pressure = 0.0d0
end subroutine parameterset
