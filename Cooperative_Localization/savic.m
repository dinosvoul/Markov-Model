function [Zazo] = savic(R,Distance)

Zazo=exp(-Distance.^2/(2*R^2));
end