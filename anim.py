import pygame
import sys
import numpy as np

# Import parameters and solvers from your math file
from freudenstein_newton_raphson import L1, L2, L3, L4, newton_raphson

# --- Pygame Setup ---
pygame.init()
WIDTH, HEIGHT = 900, 650
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("4-Bar Linkage: Multi-Slider Diagnostic Tool")
clock = pygame.time.Clock()

# Colors matching your vector diagram
BACKGROUND = (255, 255, 255)
GRID_COLOR = (240, 240, 240)
COLOR_A1 = (0, 128, 0)       # Green vector a1
COLOR_A2 = (255, 0, 0)       # Red vector a2
COLOR_A3 = (0, 0, 255)       # Blue vector a3
COLOR_A4 = (235, 50, 235)    # Magenta vector a4
TEXT_COLOR = (15, 23, 42)

# UI Colors
SLIDER_BG = (220, 225, 230)
FILL_BETA = (235, 50, 235)    # Magenta theme for beta
FILL_X0 = (100, 116, 139)     # Slate theme for x0 guess

# Fonts
font = pygame.font.SysFont("Arial", 16, bold=True)
small_font = pygame.font.SysFont("Arial", 14)

# Scale and Base Origin for the Left Ground Joint (1 cm = 25 pixels)
SCALE = 25.0
LEFT_GROUND = (WIDTH // 2 - 150, HEIGHT // 2 + 50)

# Slider Configurations (0 to 2*pi range)
BETA_MIN, BETA_MAX = 0.0, 2 * np.pi
X0_MIN, X0_MAX = 0.0, 2 * np.pi

current_beta = np.pi / 3        # Default starting Beta (60 degrees)
current_x0 = np.pi / 6          # Default starting x0 guess (30 degrees)

# Slider Positions
SLIDER_X = 150
SLIDER_W = 600
SLIDER_H = 8
HANDLE_R = 12

SLIDER_BETA_Y = 510
SLIDER_X0_Y = 570

dragging_beta = False
dragging_x0 = False

def get_joints_from_math(beta_val, x0_val):
    # 1. Left Ground Joint (Origin of x-y axes)
    Lx, Ly = LEFT_GROUND
    
    # 2. Right Ground Joint (End of vector a1, shifted horizontally by L1)
    Rx = Lx + L1 * SCALE
    Ry = Ly
    
    # 3. Top-Right Joint (Extending relative to Right Ground via beta)
    TRx = Rx + L4 * SCALE * np.cos(beta_val)
    TRy = Ry - L4 * SCALE * np.sin(beta_val)  # Inverted Y for screen graphics
    
    # Run your exact Newton-Raphson function using the manually selected x0 seed
    x_sol, iterations, converged = newton_raphson(x0_val, beta_val)
    
    if converged:
        alpha_angle = x_sol
        # 4. Top-Left Joint (Calculated from Left Ground using solved alpha angle)
        TLx = Lx + L2 * SCALE * np.cos(alpha_angle)
        TLy = Ly - L2 * SCALE * np.sin(alpha_angle)
    else:
        TLx, TLy = None, None
        alpha_angle = x_sol

    return (Lx, Ly), (Rx, Ry), (TLx, TLy), (TRx, TRy), alpha_angle, iterations, converged


# --- Main Loop ---
running = True
while running:
    mouse_pos = pygame.mouse.get_pos()
    
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
            
        elif event.type == pygame.MOUSEBUTTONDOWN:
            # Handle Beta Slider Check
            handle_beta_x = SLIDER_X + int((current_beta - BETA_MIN) / (BETA_MAX - BETA_MIN) * SLIDER_W)
            if (mouse_pos[0] - handle_beta_x)**2 + (mouse_pos[1] - SLIDER_BETA_Y)**2 <= HANDLE_R**2:
                dragging_beta = True
                
            # Handle x0 Slider Check
            handle_x0_x = SLIDER_X + int((current_x0 - X0_MIN) / (X0_MAX - X0_MIN) * SLIDER_W)
            if (mouse_pos[0] - handle_x0_x)**2 + (mouse_pos[1] - SLIDER_X0_Y)**2 <= HANDLE_R**2:
                dragging_x0 = True
                
        elif event.type == pygame.MOUSEBUTTONUP:
            dragging_beta = False
            dragging_x0 = False

    # Update dynamic slider values upon dragging
    if dragging_beta:
        relative_x = max(0, min(mouse_pos[0] - SLIDER_X, SLIDER_W))
        current_beta = BETA_MIN + (relative_x / SLIDER_W) * (BETA_MAX - BETA_MIN)
        
    if dragging_x0:
        relative_x = max(0, min(mouse_pos[0] - SLIDER_X, SLIDER_W))
        current_x0 = X0_MIN + (relative_x / SLIDER_W) * (X0_MAX - X0_MIN)

    # Compute joint coordinates using your analytical conditions
    j_LeftG, j_RightG, j_TopL, j_TopR, alpha_solved, iters, converged = get_joints_from_math(current_beta, current_x0)

    # --- Render Stage ---
    screen.fill(BACKGROUND)
    
    # Structural Gridlines
    for x in range(0, WIDTH, 40): pygame.draw.line(screen, GRID_COLOR, (x, 0), (x, HEIGHT), 1)
    for y in range(0, HEIGHT, 40): pygame.draw.line(screen, GRID_COLOR, (0, y), (WIDTH, y), 1)

    # Always Draw Ground Vector a1
    pygame.draw.line(screen, COLOR_A1, j_LeftG, j_RightG, 5)
    pygame.draw.circle(screen, TEXT_COLOR, j_LeftG, 6, 2)
    pygame.draw.circle(screen, TEXT_COLOR, j_RightG, 6, 2)

    # Always Draw Right Side Rocker Vector a4 (Controlled explicitly by beta slider)
    pygame.draw.line(screen, COLOR_A4, j_RightG, (int(j_TopR[0]), int(j_TopR[1])), 5)
    pygame.draw.circle(screen, TEXT_COLOR, (int(j_TopR[0]), int(j_TopR[1])), 6, 2)

    # ONLY draw the left side parameters if the custom guess allowed NR to converge cleanly
    if converged and j_TopL is not None:
        # Draw Crank Vector a2 (Using solved alpha angle)
        pygame.draw.line(screen, COLOR_A2, j_LeftG, j_TopL, 5)
        # Draw Coupler Vector a3 (Connecting Top-Left to Top-Right rigidly)
        pygame.draw.line(screen, COLOR_A3, j_TopL, j_TopR, 5)
        
        # Draw remaining joint dot
        pygame.draw.circle(screen, TEXT_COLOR, (int(j_TopL[0]), int(j_TopL[1])), 6, 2)
        
        # Structural Vector Labels
        screen.blit(font.render("a2", True, COLOR_A2), ((j_LeftG[0] + j_TopL[0])//2 - 20, (j_LeftG[1] + j_TopL[1])//2))
        screen.blit(font.render("a3", True, COLOR_A3), ((j_TopL[0] + j_TopR[0])//2, (j_TopL[1] + j_TopR[1])//2 - 20))

    # Base Vector Labels
    screen.blit(font.render("a1", True, COLOR_A1), ((j_LeftG[0] + j_RightG[0])//2, j_LeftG[1] + 15))
    screen.blit(font.render("a4", True, COLOR_A4), ((j_TopR[0] + j_RightG[0])//2 + 15, (j_TopR[1] + j_RightG[1])//2))


    # --- Render Slider 1: Beta (0° - 360°) ---
    pygame.draw.rect(screen, SLIDER_BG, (SLIDER_X, SLIDER_BETA_Y, SLIDER_W, SLIDER_H), border_radius=4)
    fill_beta_w = int((current_beta - BETA_MIN) / (BETA_MAX - BETA_MIN) * SLIDER_W)
    pygame.draw.rect(screen, FILL_BETA, (SLIDER_X, SLIDER_BETA_Y, fill_beta_w, SLIDER_H), border_radius=4)
    
    h_beta_x = SLIDER_X + fill_beta_w
    pygame.draw.circle(screen, FILL_BETA, (h_beta_x, SLIDER_BETA_Y), HANDLE_R)
    pygame.draw.circle(screen, (255, 255, 255), (h_beta_x, SLIDER_BETA_Y), 4)
    
    screen.blit(small_font.render("Beta (β)", True, TEXT_COLOR), (40, SLIDER_BETA_Y - 5))
    screen.blit(small_font.render("0°", True, TEXT_COLOR), (SLIDER_X, SLIDER_BETA_Y + 12))
    screen.blit(small_font.render("360°", True, TEXT_COLOR), (SLIDER_X + SLIDER_W - 30, SLIDER_BETA_Y + 12))


    # --- Render Slider 2: Initial Guess x0 (0° - 360°) ---
    pygame.draw.rect(screen, SLIDER_BG, (SLIDER_X, SLIDER_X0_Y, SLIDER_W, SLIDER_H), border_radius=4)
    fill_x0_w = int((current_x0 - X0_MIN) / (X0_MAX - X0_MIN) * SLIDER_W)
    pygame.draw.rect(screen, FILL_X0, (SLIDER_X, SLIDER_X0_Y, fill_x0_w, SLIDER_H), border_radius=4)
    
    h_x0_x = SLIDER_X + fill_x0_w
    pygame.draw.circle(screen, FILL_X0, (h_x0_x, SLIDER_X0_Y), HANDLE_R)
    pygame.draw.circle(screen, (255, 255, 255), (h_x0_x, SLIDER_X0_Y), 4)
    
    screen.blit(small_font.render("Guess (x₀)", True, TEXT_COLOR), (40, SLIDER_X0_Y - 5))
    screen.blit(small_font.render("0°", True, TEXT_COLOR), (SLIDER_X, SLIDER_X0_Y + 12))
    screen.blit(small_font.render("360°", True, TEXT_COLOR), (SLIDER_X + SLIDER_W - 30, SLIDER_X0_Y + 12))


    # --- Information Dashboard ---
    beta_deg = np.degrees(current_beta)
    x0_deg = np.degrees(current_x0)
    alpha_deg = np.degrees(alpha_solved)
    
    screen.blit(font.render(f"Input Beta (β): {beta_deg:.1f}°", True, TEXT_COLOR), (40, 30))
    screen.blit(font.render(f"Initial Seed (x₀): {x0_deg:.1f}°", True, TEXT_COLOR), (40, 55))
    screen.blit(font.render(f"NR Output Alpha (α): {alpha_deg:.1f}°", True, TEXT_COLOR), (40, 80))
    screen.blit(font.render(f"Iterations Taken: {iters}", True, TEXT_COLOR), (40, 105))
    
    status_str = "CONVERGED" if converged else "FAILED TO CONVERGE (Render Halted)"
    status_color = (16, 185, 129) if converged else (220, 38, 38)
    screen.blit(font.render(f"Solver Status: {status_str}", True, status_color), (40, 130))

    pygame.display.flip()
    clock.tick(60)

pygame.quit()
sys.exit()