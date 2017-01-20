// Copyright (C) 2017 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup spop_trimat
//! @{



template<typename T1>
inline
void
spop_trimat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_trimat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename   T1::elem_type  eT;
  typedef typename umat::elem_type ueT;
  
  const SpProxy<T1> p(in.m);
  
  arma_debug_check( (p.get_n_rows() != p.get_n_cols()), "trimatu()/trimatl(): given matrix must be square sized" );
  
  const bool upper = (in.aux_uword_a == 0);
  
  const uword N = p.get_n_nonzero();
  
  if(N == uword(0))
    {
    out.zeros(p.get_n_rows(), p.get_n_cols());
    return;
    }
  
  umat out_locs(2, N);
  
  Col<eT> out_vals(N);
  eT*     out_vals_ptr = out_vals.memptr();
  
  uword out_count = 0;
  
  typename SpProxy<T1>::const_iterator_type it = p.begin();
  
  if(upper)
    {
    // upper triangular: copy the diagonal and the elements above the diagonal
    
    for(uword in_count = 0; in_count < N; ++in_count)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      if(row <= col)
        {
        ueT* out_locs_ptr = out_locs.colptr(out_count);
        
        out_locs_ptr[0] = row;
        out_locs_ptr[1] = col;
        
        out_vals_ptr[out_count] = (*it);
        
        out_count++;
        }
      
      ++it;
      }
    }
  else
    {
    // lower triangular: copy the diagonal and the elements below the diagonal
    
    for(uword in_count = 0; in_count < N; ++in_count)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      if(row >= col)
        {
        ueT* out_locs_ptr = out_locs.colptr(out_count);
        
        out_locs_ptr[0] = row;
        out_locs_ptr[1] = col;
        
        out_vals_ptr[out_count] = (*it);
        
        out_count++;
        }
      
      ++it;
      }
    }
  
  SpMat<eT> tmp(out_locs.head_cols(out_count), out_vals.head(out_count), p.get_n_rows(), p.get_n_cols());
  
  out.steal_mem(tmp);
  }



//! @}
