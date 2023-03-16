/****** pqpqp_force_gradient.c  -- ****************/
/* MIMD version 7 */
/* Define functions that are useful for PQPQP and QPQPQ and
   force-gradient integration                              */
/* A.W-L. - adding integrators 08/2021 
 */

#include "ks_imp_includes.h"    /* definitions files and prototypes */
#ifdef HAVE_QUDA
#include "../include/generic_quda.h"
#endif

// These inner QPQPQ, PQPQP and FGI functions only update the gauge fields
// They are for doing multi-steps for the gauge field per step for the fermions

void update_inner_qpqpq( Real tau, int steps, Real lambda, int q_inner) {

    Real dtau = tau / steps;

    if(q_inner == 0){
        /* regular Q update, no need for inner_steps or lambda */
        update_u      ( tau );
    } else if(q_inner == 1){
        /* do a QPQ update in a loop of inner steps */
        for(int step=1; step <= steps; step++){
            update_u      ( dtau *0.5 );
            update_h_gauge( dtau      );
            update_u      ( dtau *0.5 );
        }
    } else if(q_inner == 2){
        for(int step=1; step <= steps; step++){
            /* update U's and H's - see header comment */
            update_u      ( dtau *lambda );
            update_h_gauge( dtau *0.5 );
            update_u      ( dtau *(1-2.*lambda) );
            update_h_gauge( dtau *0.5 );
            update_u      ( dtau *lambda );
        }
    } else {
        node0_printf("The q_inner must be 0 (Q), 1 (QPQ) or 2 (QPQPQ)\n");
        terminate(1);
    }
}

void update_inner_pqpqp( Real tau, int steps, Real lambda, int q_inner) {
    Real dtau = tau / steps;

    if(q_inner == 0){
        /* regular P_Gauge then Q update, no need for inner_steps or lambda */
        update_u      (tau);
    } else if(q_inner == 1){
        /* do a PQP update in a loop of inner steps */
        for(int step=1; step <= steps; step++){
            if(step == 1){//only do first step the first iteration through loop
                update_h_gauge( 0.5 * dtau );
            }
            update_u          (       dtau );
            // double the last step to make up for the first, except for the last iteration
            if(step == steps){
                update_h_gauge( 0.5 * dtau );
            } else{
                update_h_gauge(       dtau );
            }
        }
    } else if(q_inner == 2){
        for(int step=1; step <= steps; step++){
            //only do first step the first iteration through loop
            if(step == 1){
                update_h_gauge( dtau *lambda );
            }
            update_u          ( dtau *0.5 );
            update_h_gauge    ( dtau *(1-2.*lambda) );
            update_u          ( dtau *0.5 );
            // double the last step to make up for the first, except for the last iteration
            if(step == steps){
                update_h_gauge  ( dtau *lambda );
            } else{
                update_h_gauge  ( dtau *2.0*lambda );
            }
        }
    } else {
        node0_printf("The q_inner must be 0 (Q), 1 (PQP) or 2 (PQPQP)\n");
        terminate(1);
    }
}

void update_inner_fgi(Real tau, int steps){
  
#ifdef FN
  node0_printf("update_inner_fgi: invalidate_fermion_links\n");
  invalidate_fermion_links(fn_links);
#endif
  
  Real lambda = 1.0/6.0;
  Real xi = 1.0/72.0;
  
  Real dtau = tau / steps;
  Real lambda_dt = dtau*lambda;
  Real dtauby2 = dtau / 2.0;
  Real one_minus_2lambda_dt = (1-2*lambda)*dtau;
  Real two_lambda_dt = lambda_dt*2;
  Real xi_dtdt = 2*dtau*dtau*dtau*xi;
  
  for(int step=1; step <= steps; step+=1){
    
    if(step == 1) update_h_gauge(lambda_dt);
    
    update_u(dtauby2);
    force_gradient(one_minus_2lambda_dt, xi_dtdt, NULL, 0);
    update_u(dtauby2);
    
    if(step == steps) update_h_gauge(lambda_dt);
    else              update_h_gauge(two_lambda_dt);      
    
  }
}

// Force Gradient Integrator steps
/* see Yin and Mawhinney, arXiv:1111.5059
       Kennedy, Clark, Silva, arXiv:0910.2950
       Kennedy, Silva, Clark, arXiv:1210.6600
*/

/* We want to run a PQPQP integrator 

   exp(lam dt S) exp(0.5 dt T) exp( (1-2 lam) dt S) exp(0.5 dt T) exp( lam dt S)

   which has a Shadow Hamiltonian
   
   H = T+S
   H_shadow = H + dt**2 ( (6 lam**2 -6 lam +1)/12 {S, {S, T}}
                          +(1-6 lam)/24 {T, {S, T}})
                 + dt**4

   If we choose lam = 1/6, the coefficient of the {T, {S, T}} term becomes zero

   H_shadow = H + dt**2 ( 1/72 {S, {S, T}} )
                + dt**4

   Kate et al noticed {S, {S,T}} is independent of momentum, and depends on the 2nd
   derivative of the "potential", S

   hat {S, {S, T}} = -2 S' S'' d/dp

   The Force-Gradient Integrator removes this term as well with, with lam = 1/6

   exp(dt S / 6) exp(0.5 dt T) 
      exp( 2/3 dt S - dt**3 hat{S,{S,T}} / 72 ) 
   exp(0.5 dt T) exp( dt S / 6)

   Yin and Mawhinney 

   To accomplish the Force-Gradient step
   0. make a copy of the gauge and momentum fields, and zero out p
      p_tmp = p
      u_tmp = U
      p     = 0

   1. calculate the force  ( e( dt S )
      Fj = ej(S[U])
        update_h_gauge (1.0)  - this gets force

   2. update the gauge field
      U' = exp( -dt**2 / 24 * Fj Tj) U
        update_u( - dt**2 / 24) - this uses the previous force
 

   3. calculate the force with the updated gauge field
      Fi = ei(S[U'])

   4. add force to the momentum field
      p = p + F ?

   5. restore the gauge field
      U -> U


   The structure of the integrator will be

   P  QPQPQ_inner  P_FG  QPQPQ_inner  P

   with the QPQPQ_inner functions doing multiple gauge-only updates
   and the P_FG the force_gradient step
 */

// DMH Add an arg to specify gauge action, fermion action, or rhmc action.
// action = 0 -> gauge action
// action = 1 -> ferm action
// action = 2 -> rhmc action
int force_gradient(Real eps_t, Real eps_ttt, su3_vector **multi_x, int action){
    int iters = 0;

#ifdef FN
    node0_printf("force_gradient: invalidate_fermion_links\n");
    invalidate_fermion_links(fn_links);
#endif

    // Sanity check
    if((action == 1 || action == 2) && !multi_x) {
      node0_printf("For action 1 (ferm) or 2 (rhmc) you must pass a valid `su3_vector **multi_x` array\n");
      terminate(1);
    }
    
#ifndef HAVE_QUDA // CPU field copies
    // allocate memory for gauge field copy
    su3_matrix *linkcopyXUP, *linkcopyYUP, *linkcopyZUP, *linkcopyTUP;
    linkcopyXUP = malloc(sizeof(su3_matrix)*sites_on_node);
    linkcopyYUP = malloc(sizeof(su3_matrix)*sites_on_node);
    linkcopyZUP = malloc(sizeof(su3_matrix)*sites_on_node);
    linkcopyTUP = malloc(sizeof(su3_matrix)*sites_on_node);
    // and momentum copy
    anti_hermitmat *momentumcopy;
    momentumcopy = malloc(sizeof(anti_hermitmat)*sites_on_node*4);

    // make a copy of the gauge field
    copy_gauge_field(linkcopyXUP, linkcopyYUP, linkcopyZUP, linkcopyTUP);

    //Make a copy of the momentum field and then zero it out
    copy_momentum(momentumcopy);
    zero_momentum();
#else // have GPU copy the gauge and momentum fields
    // make a copy and zero the momentum on GPU
    node0_printf("zero device momentum after making a copy\n");
    qudaMomZero();
    // make a copy of the gauge field on GPU
    node0_printf("MILC FGI QUDA: make a copy of gauge field on device\n");
    qudaGaugeCopy();
#endif

    // ultimately, we shift by p * dt * (1-2lam)
    // so, the dt**3 shift, since it is added to p, should be normalized by
    // 1 / (dt (1-2lam))
    // eps_t = dt * (1-2lam)
    // eps_ttt = dt**3 /72 * 2
    // we need to compute the gauge / ferm contributions to the momentum
    if(action == 0)               update_h_gauge(   eps_ttt / eps_t);
    else if(action == 1) iters += update_h_fermion( eps_ttt / eps_t, multi_x);
    else if(action == 2) iters += update_h_rhmc(    eps_ttt / eps_t, multi_x);
    else {
      node0_printf("The action must be 0 (gauge), 1 (ferm) or 2 (rhmc)\n");
      terminate(1);
    }

    // given the momentum kick (force), update the links to U'
    update_u( 1.0 );

    // restore the momentum so we can add our kick to it
#ifndef HAVE_QUDA
    restore_momentum(momentumcopy);
#else
    node0_printf("MILC FGI QUDA: restore device momentum\n");
    qudaMomRestore();
#endif

    // add our kick to the momentum
    // eps_t = dt * (1-2lam)
    // DMH because we restore the gauge field, we can write a thinner
    // version of this step that does not update the gauge field, rather
    // just the momenta and fermions. However, given the speed with which
    // gauge updates are performed, it's likely to not give any significant
    // speed-up.
    if(action == 0)               update_h_gauge(   eps_t);
    else if(action == 1) iters += update_h_fermion( eps_t, multi_x);
    else if(action == 2) iters += update_h_rhmc(    eps_t, multi_x);
    else {
      node0_printf("The action must be 0 (gauge), 1 (ferm) or 2 (rhmc)\n");
      terminate(1);
    }

    // restore the gauge field
#ifndef HAVE_QUDA
    restore_gauge_field(linkcopyXUP, linkcopyYUP, linkcopyZUP, linkcopyTUP);

    // free the memory
    free(linkcopyXUP);
    free(linkcopyYUP);
    free(linkcopyZUP);
    free(linkcopyTUP);
    free(momentumcopy);
#else
    node0_printf("MILC FGI QUDA: restore device gauge\n");
    qudaGaugeRestore();
#endif
    return iters;
}



// FGI helper functions

void copy_gauge_field(su3_matrix *linkcopyXUP, 
                      su3_matrix *linkcopyYUP,
                      su3_matrix *linkcopyZUP, 
                      su3_matrix *linkcopyTUP){
    /* copy link field to old_link */
    //gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
    register int i;
    register site *s;
    FORALLSITES(i,s){
        su3mat_copy(&(s->link[XUP]),linkcopyXUP+i);
        su3mat_copy(&(s->link[YUP]),linkcopyYUP+i);
        su3mat_copy(&(s->link[ZUP]),linkcopyZUP+i);
        su3mat_copy(&(s->link[TUP]),linkcopyTUP+i);
  
    }
}

void restore_gauge_field(su3_matrix *linkcopyXUP, 
                         su3_matrix *linkcopyYUP,
                         su3_matrix *linkcopyZUP, 
                         su3_matrix *linkcopyTUP){
    register int i;
    register site *s;
    FORALLSITES(i,s){
        su3mat_copy(linkcopyXUP+i, &(s->link[XUP]));
        su3mat_copy(linkcopyYUP+i, &(s->link[YUP]));
        su3mat_copy(linkcopyZUP+i, &(s->link[ZUP]));
        su3mat_copy(linkcopyTUP+i, &(s->link[TUP]));
    }
}

void zero_momentum() {
    int i,dir;
    site *s;
    anti_hermitmat* momentum;
    field_offset mom_off = F_OFFSET(mom);
    FORALLUPDIR(dir){
        FORALLSITES(i,s) {
            momentum = (anti_hermitmat *)F_PT(s,mom_off);
            (momentum+dir)->m01.real=0;
            (momentum+dir)->m01.imag=0;
            (momentum+dir)->m02.real=0;
            (momentum+dir)->m02.imag=0;
            (momentum+dir)->m12.real=0;
            (momentum+dir)->m12.imag=0;
            (momentum+dir)->m00im=0;
            (momentum+dir)->m11im=0;
            (momentum+dir)->m22im=0;
        }
    }
}

void copy_momentum(anti_hermitmat *momentumcopy) {
    int i,dir;
    site *s;
    anti_hermitmat* momentum;
    field_offset mom_off = F_OFFSET(mom);
    FORALLUPDIR(dir){
        FORALLSITES(i,s) {
            momentum = (anti_hermitmat *)F_PT(s,mom_off);
            (momentumcopy+i+dir*sites_on_node)->m01.real=(momentum+dir)->m01.real;
            (momentumcopy+i+dir*sites_on_node)->m01.imag=(momentum+dir)->m01.imag;
            (momentumcopy+i+dir*sites_on_node)->m02.real=(momentum+dir)->m02.real;
            (momentumcopy+i+dir*sites_on_node)->m02.imag=(momentum+dir)->m02.imag;
            (momentumcopy+i+dir*sites_on_node)->m12.real=(momentum+dir)->m12.real;
            (momentumcopy+i+dir*sites_on_node)->m12.imag=(momentum+dir)->m12.imag;
            (momentumcopy+i+dir*sites_on_node)->m00im=(momentum+dir)->m00im;
            (momentumcopy+i+dir*sites_on_node)->m11im=(momentum+dir)->m11im;
            (momentumcopy+i+dir*sites_on_node)->m22im=(momentum+dir)->m22im;
        }
    }
}

void restore_momentum(anti_hermitmat *momentumcopy) {
    int i,dir;
    site *s;
    anti_hermitmat* momentum;
    field_offset mom_off = F_OFFSET(mom);
    FORALLUPDIR(dir){
        FORALLSITES(i,s) {
            momentum = (anti_hermitmat *)F_PT(s,mom_off);
            (momentum+dir)->m01.real=(momentumcopy+i+dir*sites_on_node)->m01.real;
            (momentum+dir)->m01.imag=(momentumcopy+i+dir*sites_on_node)->m01.imag;
            (momentum+dir)->m02.real=(momentumcopy+i+dir*sites_on_node)->m02.real;
            (momentum+dir)->m02.imag=(momentumcopy+i+dir*sites_on_node)->m02.imag;
            (momentum+dir)->m12.real=(momentumcopy+i+dir*sites_on_node)->m12.real;
            (momentum+dir)->m12.imag=(momentumcopy+i+dir*sites_on_node)->m12.imag;
            (momentum+dir)->m00im=(momentumcopy+i+dir*sites_on_node)->m00im;
            (momentum+dir)->m11im=(momentumcopy+i+dir*sites_on_node)->m11im;
            (momentum+dir)->m22im=(momentumcopy+i+dir*sites_on_node)->m22im;
        }
    }
}
