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
        self.history = []
        self.inverse_hessian_diagonal = None
        self.prev_gradient = None
        self.prev_step = None
        self.damping = 1.0
        self.bad_step_counter = 0
        self.success_counter = 0

    def reset(self):
        """
        Resets the internal quasi-Newton state.
        """

        self.inv_hessian = None
        self.history = []
        self.inverse_hessian_diagonal = None
        self.prev_gradient = None
        self.prev_step = None
        self.damping = 1.0
        self.bad_step_counter = 0
        self.success_counter = 0

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
            'backoff': False,
            'scaled': False,
            'damping': self.damping,
        }

        grad = np.array(gradient, dtype='float64', copy=True).reshape(-1)
        diag_inv = np.array(diag_inv_hessian, dtype='float64', copy=True).reshape(-1)
        self.inverse_hessian_diagonal = diag_inv.copy()

        if self.prev_gradient is None or self.prev_gradient.size != grad.size:
            self.history = []
            self.prev_gradient = None
            self.prev_step = None
            self.damping = options['damping']
            self.bad_step_counter = 0
            self.success_counter = 0
            info['reset'] = True
        else:
            self._update_inverse_hessian_model(grad, diag_inv, options, info)

        step = -self.damping * self._apply_inverse_hessian(grad, diag_inv)

        step_norm = np.linalg.norm(step)
        max_step = options['max_step_norm']
        if step_norm > max_step and step_norm > 0.0:
            step *= max_step / step_norm
            info['scaled'] = True

        self.prev_gradient = grad
        self.prev_step = step.copy()
        info['damping'] = self.damping

        return step, info

    def _update_inverse_hessian_model(self, grad, diag_inv, options, info):
        """
        Updates or resets the inverse-Hessian model using the latest gradient.
        """

        prev_grad_norm = (np.linalg.norm(self.prev_gradient)
                          if self.prev_gradient is not None else None)
        grad_norm = np.linalg.norm(grad)

        if (self.prev_gradient is None or self.prev_step is None or
                self.prev_gradient.size != grad.size):
            self.history = []
            self.prev_gradient = None
            self.prev_step = None
            self.damping = options['damping']
            self.bad_step_counter = 0
            self.success_counter = 0
            info['reset'] = True
            return

        if prev_grad_norm is not None and prev_grad_norm > 0.0:
            ratio = grad_norm / prev_grad_norm
            if ratio > options['reset_ratio']:
                self.bad_step_counter += 1
                self.success_counter = 0
                self.damping = max(options['min_damping'],
                                   self.damping * options['damping_decay'])
                info['backoff'] = True
                if self.bad_step_counter >= options['reset_persistence']:
                    self.history = []
                    self.bad_step_counter = 0
                    info['reset'] = True
                return

            self.bad_step_counter = 0
            if ratio < options['improve_ratio']:
                self.success_counter += 1
            else:
                self.success_counter = 0

            if self.success_counter >= options['success_recovery_steps']:
                self.damping = min(options['damping'],
                                   self.damping * options['damping_growth'])
                self.success_counter = 0

        if not options['use_bfgs']:
            self.history = []
            return

        svec = self.prev_step
        yvec = grad - self.prev_gradient
        sty = np.dot(svec, yvec)
        if sty <= options['curvature_tol']:
            self.bad_step_counter += 1
            self.success_counter = 0
            self.damping = max(options['min_damping'],
                               self.damping * options['damping_decay'])
            info['backoff'] = True
            if self.bad_step_counter >= options['reset_persistence']:
                self.history = []
                self.bad_step_counter = 0
                info['reset'] = True
            return

        self.bad_step_counter = 0
        rho = 1.0 / sty
        self.history.append((svec.copy(), yvec.copy(), rho))
        max_history = max(0, int(options['history_size']))
        if max_history == 0:
            self.history = []
        elif len(self.history) > max_history:
            self.history = self.history[-max_history:]
        info['bfgs_updated'] = True

    def _apply_inverse_hessian(self, grad, diag_inv):
        """
        Applies the inverse-Hessian approximation to a vector.
        """

        if len(self.history) == 0:
            return diag_inv * grad

        qvec = grad.copy()
        alphas = []

        for svec, yvec, rho in reversed(self.history):
            alpha = rho * np.dot(svec, qvec)
            alphas.append(alpha)
            qvec = qvec - alpha * yvec

        rvec = diag_inv * qvec

        for alpha, (svec, yvec, rho) in zip(reversed(alphas), self.history):
            beta = rho * np.dot(yvec, rvec)
            rvec = rvec + svec * (alpha - beta)

        return rvec
