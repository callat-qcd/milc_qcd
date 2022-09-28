/****** pqpqp_force_gradient.c  -- ****************/
/* MIMD version 7 */
/* Define functions that are useful for PQPQP and QPQPQ and
   force-gradient integration                              */
/* A.W-L. - adding integrators 08/2021 
 */

#include "ks_imp_includes.h"    /* definitions files and prototypes */


// These inner QPQPQ and PQPQP functions only update the gauge fields
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
            if(step == 1){//only do first step the first itteration through loop
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
            //only do first step the first itteration through loop
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
        node0_printf("The q_inner must be 0 (Q), 1 (QPQ) or 2 (QPQPQ)\n");
        terminate(1);
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
   H_shaddow = H + dt**2 ( (6 lam**2 -6 lam +1)/12 {S, {S, T}}
                           +(1-6 lam)/24 {T, {S, T}})
                 + dt**4

   If we chose lam = 1/6, the coefficient of the {T, {S, T}} term becomes zero

   H_shaddow = H + dt**2 ( 1/72 {S, {S, T}} )
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

int force_gradient( Real eps_t, Real eps_ttt, su3_vector **multi_x ){
    int iters;

#ifdef FN
    invalidate_fermion_links(fn_links);
#endif
    
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

    // ultimately, we shift by p * dt * (1-2lam)
    // so, the dt**3 shift, since it is added to p, should be normalized by
    // 1 / (dt (1-2lam))
    // eps_t = dt * (1-2lam)
    // eps_ttt = dt**3 /72 * 2
    // we need to compute the gauge and ferm contributions to the momentum
    iters += update_h_rhmc( eps_ttt / eps_t, multi_x);

    // given the momentum kick (force), update the links to U'
    update_u( 1.0 );

    // restore the momentum so we can add our kick to it
    restore_momentum(momentumcopy);

    // add our kick to the momentum
    // eps_t = dt * (1-2lam)
    iters += update_h_rhmc( eps_t, multi_x );

    // restore the gauge field
    restore_gauge_field(linkcopyXUP, linkcopyYUP, linkcopyZUP, linkcopyTUP);

    // free the memory
    free(linkcopyXUP);
    free(linkcopyYUP);
    free(linkcopyZUP);
    free(linkcopyTUP);
    free(momentumcopy);
    
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
