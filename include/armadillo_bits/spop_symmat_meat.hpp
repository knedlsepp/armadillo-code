// Copyright (C) 2017 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup spop_symmat
//! @{



template<typename T1>
inline
void
spop_symmat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_symmat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename   T1::elem_type  eT;
  typedef typename umat::elem_type ueT;
  
  const SpProxy<T1> p(in.m);
  
  arma_debug_check( (p.get_n_rows() != p.get_n_cols()), "symmatu()/symmatl(): given matrix must be square sized" );
  
  const bool upper = (in.aux_uword_a == 0);
  
  const uword N = p.get_n_nonzero();
  
  if(N == uword(0))
    {
    out.zeros(p.get_n_rows(), p.get_n_cols());
    return;
    }
  
  umat out_locs(2, 2*N);   // 2*N for worst case scenario
  
  Col<eT> out_vals(2*N);
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
      
      if(row < col)
        {
        ueT* out_locs_ptr_a = out_locs.colptr(out_count  );
        ueT* out_locs_ptr_b = out_locs.colptr(out_count+1);
        
        out_locs_ptr_a[0] = row;
        out_locs_ptr_a[1] = col;
        
        out_locs_ptr_b[0] = col;
        out_locs_ptr_b[1] = row;
        
        const eT val = (*it);
        
        out_vals_ptr[out_count  ] = val;
        out_vals_ptr[out_count+1] = val;
        
        out_count += 2;
        }
      else
      if(row == col)
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
      
      if(row > col)
        {
        ueT* out_locs_ptr_a = out_locs.colptr(out_count  );
        ueT* out_locs_ptr_b = out_locs.colptr(out_count+1);
        
        out_locs_ptr_a[0] = row;
        out_locs_ptr_a[1] = col;
        
        out_locs_ptr_b[0] = col;
        out_locs_ptr_b[1] = row;
        
        const eT val = (*it);
        
        out_vals_ptr[out_count  ] = val;
        out_vals_ptr[out_count+1] = val;
        
        out_count += 2;
        }
      else
      if(row == col)
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



template<typename T1>
inline
void
spop_symmat_cx::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_symmat_cx>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename   T1::elem_type  eT;
  typedef typename umat::elem_type ueT;
  
  const SpProxy<T1> p(in.m);
  
  arma_debug_check( (p.get_n_rows() != p.get_n_cols()), "symmatu()/symmatl(): given matrix must be square sized" );
  
  const bool upper   = (in.aux_uword_a == 0);
  const bool do_conj = (in.aux_uword_b == 1);
  
  const uword N = p.get_n_nonzero();
  
  if(N == uword(0))
    {
    out.zeros(p.get_n_rows(), p.get_n_cols());
    return;
    }
  
  umat out_locs(2, 2*N);   // 2*N for worst case scenario
  
  Col<eT> out_vals(2*N);
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
      
      if(row < col)
        {
        ueT* out_locs_ptr_a = out_locs.colptr(out_count  );
        ueT* out_locs_ptr_b = out_locs.colptr(out_count+1);
        
        out_locs_ptr_a[0] = row;
        out_locs_ptr_a[1] = col;
        
        out_locs_ptr_b[0] = col;
        out_locs_ptr_b[1] = row;
        
        const eT val = (*it);
        
        out_vals_ptr[out_count  ] = val;
        out_vals_ptr[out_count+1] = (do_conj) ? std::conj(val) : val;
        
        out_count += 2;
        }
      else
      if(row == col)
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
      
      if(row > col)
        {
        ueT* out_locs_ptr_a = out_locs.colptr(out_count  );
        ueT* out_locs_ptr_b = out_locs.colptr(out_count+1);
        
        out_locs_ptr_a[0] = row;
        out_locs_ptr_a[1] = col;
        
        out_locs_ptr_b[0] = col;
        out_locs_ptr_b[1] = row;
        
        const eT val = (*it);
        
        out_vals_ptr[out_count  ] = val;
        out_vals_ptr[out_count+1] = (do_conj) ? std::conj(val) : val;
        
        out_count += 2;
        }
      else
      if(row == col)
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
