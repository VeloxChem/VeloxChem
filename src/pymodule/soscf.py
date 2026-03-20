#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np


class Soscf:
    """
    Simple SOSCF accelerator based on occupied-virtual rotations.

    The initial inverse Hessian is built from orbital energy differences.
    Optionally, a dense inverse-BFGS update is accumulated in the orbital
    rotation space.
    """

    def __init__(self):
        self.inv_hessian = None
        self.prev_gradient = None
        self.prev_step = None
        self.damping = 1.0

    def reset(self):
        """
        Resets the internal quasi-Newton state.
        """

        self.inv_hessian = None
        self.prev_gradient = None
        self.prev_step = None
        self.damping = 1.0

    def compute_step(self, gradient, diag_inv_hessian, options):
        """
        Computes the SOSCF step and updates the internal quasi-Newton state.

        :param gradient:
            The orbital gradient vector.
        :param diag_inv_hessian:
            The diagonal inverse Hessian estimate.
        :param options:
            The SOSCF options dictionary.

        :return:
            The step vector and bookkeeping information.
        """

        info = {
            'bfgs_updated': False,
            'reset': False,
            'scaled': False,
            'damping': self.damping,
        }

        grad = np.array(gradient, dtype='float64', copy=True).reshape(-1)
        diag_inv = np.array(diag_inv_hessian, dtype='float64', copy=True).reshape(-1)

        if self.inv_hessian is None or self.inv_hessian.shape != (grad.size, grad.size):
            self.inv_hessian = np.diag(diag_inv)
            self.prev_gradient = None
            self.prev_step = None
            self.damping = options['damping']
            info['reset'] = True
        else:
            self._update_inverse_hessian(grad, diag_inv, options, info)

        step = -self.damping * np.dot(self.inv_hessian, grad)

        step_norm = np.linalg.norm(step)
        max_step = options['max_step_norm']
        if step_norm > max_step and step_norm > 0.0:
            step *= max_step / step_norm
            info['scaled'] = True

        self.prev_gradient = grad
        self.prev_step = step.copy()
        info['damping'] = self.damping

        return step, info

    def _update_inverse_hessian(self, grad, diag_inv, options, info):
        """
        Updates or resets the inverse Hessian using the latest gradient.
        """

        prev_grad_norm = (np.linalg.norm(self.prev_gradient)
                          if self.prev_gradient is not None else None)
        grad_norm = np.linalg.norm(grad)

        if (self.prev_gradient is None or self.prev_step is None or
                self.prev_gradient.size != grad.size):
            self.inv_hessian = np.diag(diag_inv)
            self.prev_gradient = None
            self.prev_step = None
            self.damping = options['damping']
            info['reset'] = True
            return

        if prev_grad_norm is not None and prev_grad_norm > 0.0:
            ratio = grad_norm / prev_grad_norm
            if ratio > options['reset_ratio']:
                self.inv_hessian = np.diag(diag_inv)
                self.damping = max(options['min_damping'],
                                   self.damping * options['damping_decay'])
                info['reset'] = True
                return

            if ratio < options['improve_ratio']:
                self.damping = min(options['damping'],
                                   self.damping * options['damping_growth'])

        if not options['use_bfgs']:
            self.inv_hessian = np.diag(diag_inv)
            return

        svec = self.prev_step
        yvec = grad - self.prev_gradient
        sty = np.dot(svec, yvec)
        if sty <= options['curvature_tol']:
            self.inv_hessian = np.diag(diag_inv)
            info['reset'] = True
            return

        dim = grad.size
        ident = np.eye(dim)
        rho = 1.0 / sty
        left = ident - rho * np.outer(svec, yvec)
        right = ident - rho * np.outer(yvec, svec)
        updated = np.linalg.multi_dot([left, self.inv_hessian, right])
        updated += rho * np.outer(svec, svec)

        # Keep the inverse Hessian symmetric and blend the diagonal
        # preconditioner back in to limit drift.
        updated = 0.5 * (updated + updated.T)
        updated += options['diag_blend'] * np.diag(diag_inv)

        self.inv_hessian = updated
        info['bfgs_updated'] = True
